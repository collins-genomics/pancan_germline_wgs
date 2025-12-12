# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Refine SV genotypes based on flanking SNV genotypes

# Note: this workflow assumes variant callset properties follow GATK-SV and
# GATK-HC conventions for short-read WGS. It has not been tested on other 
# callers (e.g., DeepVariant) or other technologies (e.g., long-read WGS)


version 1.0


import "Utilities.wdl" as Utils


workflow RefineSvGenotypesWithSnvs {
  input {
    File sv_vcf
    File sv_vcf_idx

    # Two ways to provide SNV VCF information: as arrays for VCF & indexes, or
    # as a two-column .tsv with URIs for VCF and index. If both are provided,
    # the array-style inputs will be used. Also note that these VCFs don't need
    # to exclusively include SNVs, but all non-SNVs will be dropped for genotyping
    Array[File]? snv_vcfs
    Array[File]? snv_vcf_idxs
    File? snv_vcf_info_tsv

    # SV filtering parameters for GT imputation training
    Float min_sv_af = 0.05               # Minimum SV AF to be regenotyped. Max is 1 - this
    Int min_sv_ac = 20                   # Minimum SV AC to be regenotyped. Max is max(AN) - this
    Int min_an = 100                     # Minimum AN to be regenotyped for both SNVs and SVs
    Boolean mask_training_sv_gts = true  # Should SV GTs be filtered for quality before model training?
    String sv_training_mask_field = "SL" # FORMAT field to use for masking low-quality genotypes prior to training

    # SNV filtering parameters
    Int breakpoint_buffer_bp = 5000      # SNVs closer than this distance to each breakpoint will not be included
    Int breakpoint_window_bp = 100000    # SNVs farther than this distance + buffer from each breakpoint will not be included
    Float min_snv_call_rate = 0.95       # Minimum call rate for SNVs to be included
    Float snv_freq_scalar = 5            # Frequency control parameter for SNVs, defined as AF ~ [min_sv_af / this, min_sv_af * this]
    File? snv_exclusion_bed              # SNVs overlapping this BED file will be excluded

    # Imputation parameters
    Int max_snps_per_flank = 10          # Max number of SNPs to include per flank
    Float min_ld_r2 = 0.1                # Minimum LD between SNV & SV to permit for imputation
    String ref_build = "hg38"            # Plink-styled reference indicator

    # Parallelization options
    Int svs_per_shard = 200
    Int snv_vcfs_per_shard = 100

    File genome_file                     # BEDTools-style .genome file

    String output_prefix

    String g2c_analysis_docker
    String tmp_dev_docker
    String linux_docker
  }

  # Determine method of SNV VCF input
  if ( defined(snv_vcfs) && defined(snv_vcf_idxs) ) {
    call WriteSnvInfo {
      input:
        vcf_uris = select_first([snv_vcfs, []]),
        tbi_uris = select_first([snv_vcf_idxs, []]),
        output_prefix = output_prefix,
        docker = linux_docker
    }
  }
  File snv_info_tsv = select_first(select_all([WriteSnvInfo.snv_info_tsv, snv_vcf_info_tsv]))
  if ( !defined(snv_vcfs) ) {
    call GetFirstSnvVcf {
      input:
        tsv = select_first(select_all([snv_vcf_info_tsv])),
        docker = linux_docker
    }
  }
  File first_snv_vcf = if defined(snv_vcfs) 
                       then select_first([snv_vcfs, []])[0] 
                       else select_first(select_all([GetFirstSnvVcf.vcf_uri]))
  File first_snv_idx = if defined(snv_vcf_idxs) 
                       then select_first([snv_vcf_idxs, []])[0] 
                       else select_first(select_all([GetFirstSnvVcf.tbi_uri]))

  # Gather list of matching samples between snv_vcfs and sv_vcf
  call Utils.StreamSamplesFromVcfHeader as GetSvSamples {
    input:
      vcf = sv_vcf,
      vcf_idx = sv_vcf_idx,
      bcftools_docker = g2c_analysis_docker
  }
  call Utils.StreamSamplesFromVcfHeader as GetSnvSamples {
    input:
      vcf = first_snv_vcf,
      vcf_idx = first_snv_idx,
      bcftools_docker = g2c_analysis_docker
  }
  Array[File] raw_sample_lists = select_all([GetSvSamples.sample_list, GetSnvSamples.sample_list])
  call Utils.IntersectTextFiles as FindSharedSamples {
    input:
      files = raw_sample_lists,
      outfile = output_prefix + ".shared_samples.list",
      docker = linux_docker
  }
  Int max_sv_ac = floor((2 * GetSvSamples.n_samples) - min_sv_ac)

  # Separate initial qualifying and non-qualifying SVs
  call SplitSvs {
    input:
      vcf = sv_vcf,
      vcf_idx = sv_vcf_idx,
      min_af = min_sv_af,
      max_af = 1 - min_sv_af,
      min_ac = min_sv_ac,
      max_ac = max_sv_ac,
      min_an = min_an,
      output_prefix = basename(sv_vcf, ".vcf.gz"),
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Shard qualifying SV VCF for parallel processing
  call Utils.ShardVcf as ShardTargetSvs {
    input:
      vcf = SplitSvs.target_sv_vcf,
      vcf_idx = SplitSvs.target_sv_vcf_idx,
      records_per_shard = svs_per_shard,
      bcftools_docker = g2c_analysis_docker,
      n_preemptible = 1
  }

  # Shard SNV VCF information
  call Utils.ShardTextFile as ShardVcfInfo {
    input:
      input_file = snv_info_tsv,
      lines_per_split = snv_vcfs_per_shard,
      out_prefix = output_prefix + ".snv_vcf_info",
      g2c_analysis_docker = g2c_analysis_docker
  }
  Array[File] vcf_info_chunks = ShardVcfInfo.shards

  # Process each qualifying SV VCF in parallel
  scatter ( vcf_info in zip(ShardTargetSvs.vcf_shards, ShardTargetSvs.vcf_shard_idxs) ) {

    String shard_prefix = basename(vcf_info.left, ".vcf.gz")

    # Apply minimal SV quality filters to enrich for true-positive SV genotypes
    # Also ensure all SVs still meet frequency criteria after filtering
    if ( mask_training_sv_gts ) {
      call MaskSvGts {
        input:
          vcf = vcf_info.left,
          vcf_idx = vcf_info.right,
          sv_training_mask_field = sv_training_mask_field,
          min_af = min_sv_af,
          max_af = 1 - min_sv_af,
          min_ac = min_sv_ac,
          max_ac = max_sv_ac,
          min_an = min_an,
          output_prefix = output_prefix + ".sl_masked",
          g2c_analysis_docker = tmp_dev_docker
      }
    }
    File sv_training_vcf = select_first([MaskSvGts.filtered_vcf, vcf_info.left])
    File sv_training_vcf_idx = select_first([MaskSvGts.filtered_vcf_idx, vcf_info.right])

    # Only proceed if >0 SVs were retained after quality masking
    Int n_sv_in_shard = select_first([MaskSvGts.n_eligible_svs, snv_vcfs_per_shard])
    if ( n_sv_in_shard > 0 ) {
      # Define SNV query intervals for this shard
    call DefineQueryIntervals {
      input:
        sv_vcf = sv_training_vcf,
        sv_vcf_idx = sv_training_vcf_idx,
        genome_file = genome_file,
        breakpoint_buffer_bp = breakpoint_buffer_bp,
        breakpoint_window_bp = breakpoint_window_bp,
        snv_freq_scalar = snv_freq_scalar,
        snv_exclusion_bed = snv_exclusion_bed,
        output_prefix = shard_prefix,
        bcftools_docker = g2c_analysis_docker
    }

    # Scatter over SNV VCF chunks and filter SNVs
    scatter ( i in range(length(vcf_info_chunks)) ) {
      call QuerySnvs {
        input:
          snv_info_tsv = vcf_info_chunks[i],
          query_intervals = DefineQueryIntervals.query_intervals,
          samples_list = FindSharedSamples.intersection_file,
          min_ac = floor(min_sv_ac / snv_freq_scalar),
          min_an = min_an,
          max_ncr = 1 - min_snv_call_rate,
          output_prefix = shard_prefix + ".chunk_" + i
      }
    }
    call Utils.ConcatVcfs as ConcatSnvs {
      input:
        vcfs = QuerySnvs.snv_vcf,
        vcf_idxs = QuerySnvs.snv_vcf_idx,
        out_prefix = shard_prefix + ".eligible_snvs",
        bcftools_concat_options = "--allow-overlaps --remove-duplicates",
        bcftools_docker = g2c_analysis_docker
    }

    # Compute LD for each SV, extract AD matrixes, fit regression model, and predict GTs for all samples
    # TODO: need to update this to train on filtered SV VCF but apply to full cohort of SNV GTs
    call ImputeSvs {
      input:
        sv_vcf = sv_training_vcf,
        sv_vcf_idx = sv_training_vcf_idx,
        snv_vcf = ConcatSnvs.merged_vcf,
        snv_vcf_idx = ConcatSnvs.merged_vcf_idx,
        training_samples_list = FindSharedSamples.intersection_file,
        breakpoint_buffer_bp = breakpoint_buffer_bp,
        breakpoint_window_bp = breakpoint_window_bp,
        snv_freq_scalar = snv_freq_scalar,
        min_ld_r2 = min_ld_r2,
        ref_build = ref_build,
        max_snps_per_flank = max_snps_per_flank,
        output_prefix = shard_prefix,
        g2c_analysis_docker = tmp_dev_docker
    }

    # Update SV GTs
    # TODO: implement this
    # NOTE: must apply to vcf_info.left, *not* quality-filtered VCF
    # - If SNV-predicted GQ > GATK-SV GQ, return (sample, SV ID, GT, GQ) to be updated in SV VCF
    }
    # TODO: add select_first() statement here to retain updated SV VCF (if generated), or initial SV VCF (if not)
  }

  # Concatenate all updated SV VCFs with the passthrough VCF
  # TODO: implement this

  output {}
}


# Defines interval information for SNV querying
task DefineQueryIntervals {
  input {
    File sv_vcf
    File sv_vcf_idx
    File genome_file

    Int breakpoint_buffer_bp
    Int breakpoint_window_bp
    Float snv_freq_scalar
    File? snv_exclusion_bed

    String output_prefix

    String bcftools_docker
  }

  String outfile = output_prefix + ".snv_query_intervals.tsv"

  String excl_cmd = if defined(snv_exclusion_bed) 
                    then "| bedtools subtract -a - -b " + basename(select_first([snv_exclusion_bed, ""]))
                    else ""

  Int disk_gb = ceil(2 * size(sv_vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    # Relocate SV tabix index, if necessary
    if [ "~{sv_vcf_idx}" != "~{sv_vcf}.tbi" ]; then
      ln -s ~{sv_vcf_idx} ~{sv_vcf}.tbi
    fi

    # Relocate exclusion bed to pwd, if provided
    if ~{defined(snv_exclusion_bed)}; then
      ln -s ~{default="" snv_exclusion_bed} .
    fi

    # To ensure GATK strict interval compliance, we trim all intervals to the bounds of all chromosomes
    awk -v OFS="\t" '{ print $1, $2-1, $2+100000000 }' ~{genome_file} > cliff.bed

    # Build query intervals
    bcftools query \
      -f '%CHROM\t%POS\t%INFO/END\t%INFO/AF\n' \
      ~{sv_vcf} \
    | awk -v FS="\t" -v OFS="\t" -v scalar=~{snv_freq_scalar} \
        -v buffer=~{breakpoint_buffer_bp} -v window=~{breakpoint_window_bp} \
        '{ print $1, $2-buffer-window, $2-buffer, $4/scalar, $4*scalar"\n"\
                 $1, $3+buffer, $3+buffer+window, $4/scalar, $4*scalar }' \
    | awk -v FS="\t" -v OFS="\t" '{ if ($2<0) $2=0; if ($5>1) $5=1; print }' \
    ~{excl_cmd} \
    | bedtools subtract -a - -b cliff.bed \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | bedtools merge -i - -d 5000 -c 4,5 -o min,max \
    | awk -v OFS="\t" '{ print NR, $0 }' \
    > ~{outfile}
  >>>

  output {
    File query_intervals = outfile
  }

  runtime {
    docker: bcftools_docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk " + disk_gb +" HDD"
    preemptible: 1
    maxRetries: 1
  }
}


# Extracts the first SNV VCF URI from a .tsv of VCF info
task GetFirstSnvVcf {
  input {
    File tsv
    String docker
  }

  command <<<
    set -eu -o pipefail

    sed -n '1p' ~{tsv} > first_info.tsv
    awk -v FS="\t" '{ print $1 }' first_info.tsv > vcf.list
    awk -v FS="\t" '{ print $2 }' first_info.tsv > tbi.list
  >>>

  output {
    String vcf_uri = read_string("vcf.list")
    String tbi_uri = read_string("tbi.list")
  }

  runtime {
    docker: docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk 25 HDD"
    preemptible: 1
    maxRetries: 1
  }
}


# Combines curated SNVs & SVs to impute SV genotypes for all eligible samples
task ImputeSvs {
  input {
    File sv_vcf
    File sv_vcf_idx
    
    File snv_vcf
    File snv_vcf_idx
    
    File training_samples_list
    Int breakpoint_buffer_bp
    Int breakpoint_window_bp
    Float snv_freq_scalar

    Float min_ld_r2
    String ref_build
    Int max_snps_per_flank
    
    String output_prefix

    String g2c_analysis_docker

    Float mem_gb = 15.5
    Int n_cpu = 8
    Int n_preemptible = 1
  }

  Int disk_gb = ceil(2 * size([sv_vcf, snv_vcf], "GB")) + 10

  Int bcftools_threads = (2 * n_cpu) - 1

  command <<<
    set -eu -o pipefail

    # Make list of SVs to process
    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/AF\n' ~{sv_vcf} \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
    > svs.bed

    # Make dummy .fam file
    awk -v OFS="\t" '{ print $1, $1, 0, 0, 0, 0 }' ~{training_samples_list} > samples.fam

    # Process each SV in serial
    while read chrom start end svid svaf; do

      echo -e "\nNow processing $svid..."

      mkdir $svid

      # Compile information for SNV filtering
      echo -e "$chrom\t$start" \
      | awk -v FS="\t" -v OFS="\t" \
        -v buffer=~{breakpoint_buffer_bp} -v window=~{breakpoint_window_bp} \
        '{ print $1, $2-buffer-window, $2-buffer }' \
      | awk -v FS="\t" -v OFS="\t" '{ if ($2<0) $2=0; print }' \
      > $svid/left.window.bed
      echo -e "$chrom\t$end" \
      | awk -v FS="\t" -v OFS="\t" \
        -v buffer=~{breakpoint_buffer_bp} -v window=~{breakpoint_window_bp} \
        '{ print $1, $2+buffer, $2+buffer+window }' \
      | awk -v FS="\t" -v OFS="\t" '{ if ($2<0) $2=0; print }' \
      > $svid/right.window.bed
      min_snv_af=$( echo $svaf | awk -v scalar=5 '{ af=($1 / scalar); if (af<0) af=0; print af }' )
      max_snv_af=$( echo $svaf | awk -v scalar=5 '{ af=($1 * scalar); if (af>1) af=1; print af }' )

      # Make VCF sandwich of (left flanking SNVs) + SV + (right flanking SNVs)
      for flank in left right; do
        bcftools view \
          --samples-file ~{training_samples_list} \
          --force-samples \
          --regions-file $svid/$flank.window.bed \
          --min-af $min_snv_af \
          --max-af $max_snv_af \
          ~{snv_vcf} \
        | /opt/pancan_germline_wgs/scripts/qc/vcf_qc/set_g2c_qc_variant_ids.py \
          --vcf-out $svid/$flank.snvs.vcf.gz
        tabix -p vcf -f $svid/$flank.snvs.vcf.gz
      done
      bcftools view \
        --samples-file ~{training_samples_list} \
        --force-samples \
        --include "ID = \"$svid\"" \
        -Oz -o $svid/sv.vcf.gz \
        ~{sv_vcf}
      tabix -p vcf -f $svid/sv.vcf.gz
      bcftools concat \
        --threads ~{bcftools_threads} \
        $svid/left.snvs.vcf.gz \
        $svid/sv.vcf.gz \
        $svid/right.snvs.vcf.gz \
      | bcftools annotate -x ^FORMAT/GT,FORMAT/AD,FORMAT/DP \
        -Oz -o $svid/sandwich.vcf.gz
      tabix -p vcf -f $svid/sandwich.vcf.gz

      # Compute all LD vs. target SV with plink
      svlen=$(( end - start ))
      if [ $svlen -lt 1 ]; then svlen=1; fi
      plink2 \
        --r2-unphased 'yes-really' 'ref-based' \
        --ld-window-kb $(( ~{breakpoint_buffer_bp} + ~{breakpoint_window_bp} + $svlen + 1 )) \
        --ld-window-r2 ~{min_ld_r2} \
        --split-par "~{ref_build}" \
        --polyploid-mode missing \
        --fam samples.fam \
        --ld-snp "$svid" \
        --vcf $svid/sandwich.vcf.gz \
        --out $svid/ld
      cut -f6,7 $svid/ld.vcor | sed '1d' | sort -nrk2,2 -k1,1V > $svid/ld.vcor.slim

      # Rank-order SNVs by LD R2 per flank and take up to N best tag SNVs from each flank
      for flank in left right; do
        bcftools query -f '%ID\n' $svid/$flank.snvs.vcf.gz > $svid/$flank.vids.list
        fgrep \
          -wf $svid/$flank.vids.list \
          $svid/ld.vcor.slim \
         | head -n ~{max_snps_per_flank} \
         | cut -f1 \
         > $svid/$flank.keep_vids.list || true
      done

      # Extract allele dosage for the SV and each flanking SNP
      cat \
        $svid/left.keep_vids.list \
        $svid/right.keep_vids.list \
        <( echo -e "$svid" ) \
      | sort -V | uniq \
      > $svid/all.keep_vids.list
      /opt/pancan_germline_wgs/scripts/variant_filtering/extract_ad_matrix.py \
        -i $svid/sandwich.vcf.gz \
        -v $svid/all.keep_vids.list \
        -o $svid/$svid.ad.tsv.gz

      # Temporary task to deloc files for development
      cp $svid/$svid.ad.tsv.gz ./

      # TODO: finish implementing this
      # - Load tag SNP AD and SV GTs into R
      # - Fit linear regression of SV AC ~ tag SNP ADs using 10-fold CV
      #   - Need to think about how to handle/prespecify train/test split (by cohort etc)
      # - Predict SV AC from tag SNPs using best-fit regression model
      # - Compute SNV-based GQ by ratio of linear distances between integer AC states (or maybe multivariate gaussian)

      # Clean up
      rm -rf $svid

    done < svs.bed
  >>>

  output {
    Array[File?] ads_tmp = glob("*.ad.tsv.gz")
  }

  runtime {
    docker: g2c_analysis_docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: n_preemptible
    maxRetries: 1
  }
}


# Apply a quality filter to all SV genotypes and retain only SVs still qualified by frequency
task MaskSvGts {
  input {
    File vcf
    File vcf_idx

    String sv_training_mask_field = "SL"
    Float min_af
    Float max_af
    Int min_ac
    Int max_ac
    Int min_an

    String output_prefix

    String g2c_analysis_docker
  }

  Int disk_gb = ceil(3.5 * size(vcf, "GB")) + 10

  String vcf_out = output_prefix + ".vcf.gz"

  command <<<
    set -eu -o pipefail

    # Apply quality mask, update AC/AN/AF, and filter again by frequency to ensure eligibility
    /opt/pancan_germline_wgs/scripts/variant_filtering/mask_sv_gts_for_regenotyping.py \
      --input-vcf ~{vcf} \
      --quality-field "~{sv_training_mask_field}" \
    | bcftools +fill-tags -- -t AC,AN,AF \
    | /opt/pancan_germline_wgs/scripts/variant_filtering/bifurcate_svs_for_regenotyping.py \
      -e "~{vcf_out}" \
      --min-af ~{min_af} \
      --max-af ~{max_af} \
      --min-ac ~{min_ac} \
      --max-ac ~{max_ac} \
      --min-an ~{min_an}
    tabix -p vcf "~{vcf_out}"

    # Count number of records after filtering to catch empty shards
    bcftools query -f '%CHROM\n' "~{vcf_out}" | wc -l > count.txt
  >>>

  output {
    File filtered_vcf = vcf_out
    File filtered_vcf_idx = "~{vcf_out}.tbi"
    Int n_eligible_svs = read_int("count.txt")
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "3.75 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 1
    maxRetries: 1
  }
}


# Extracts SNVs qualified for regenotyping from a list of SNV VCFs
task QuerySnvs {
  input {
    File snv_info_tsv
    File query_intervals
    File samples_list
    
    Int min_ac
    Int min_an
    Float max_ncr

    Int max_stream_attempts = 5

    String output_prefix

    Int disk_gb = 275
    Float mem_gb = 7.5
    Int n_cpu = 2
    Int n_preemptible = 1
    String gatk_docker = "broadinstitute/gatk:4.6.2.0"
  }

  String out_vcf = output_prefix + ".snvs.vcf.gz"
  String out_tbi = out_vcf + ".tbi"

  Int concat_threads = 2 * n_cpu
  Int sort_mem_mb = floor(1000 * (mem_gb / 3))

  command <<<
    set -eu -o pipefail

    # Get absolute minimum & maximum frequencies to permit in any interval
    global_min_af=$( cut -f5 ~{query_intervals} | sort -nk1,1 | sed -n '1p' )
    global_max_af=$( cut -f6 ~{query_intervals} | sort -nrk1,1 | sed -n '1p' )
    cut -f2-4 ~{query_intervals} > snv_query_intervals.bed

    # Make local mini-VCFs across all query intervals
    mkdir local_vcfs
    while read idx vcf_uri tbi_uri; do

      echo -e "\nQuerying $( basename $vcf_uri )..."
    
      gsutil cp "$tbi_uri" ./
      
      # Refresh token
      export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
      export GCS_TOKEN=$(gcloud auth application-default print-refresh-token 2>/dev/null || echo "")
      
      # Slice and filter qualifying SNPs
      gatk SelectVariants \
        -L snv_query_intervals.bed \
        --select-type-to-include SNP \
        --restrict-alleles-to BIALLELIC \
        --create-output-variant-index false \
        --create-output-variant-md5 false \
        -V $vcf_uri \
        --output /dev/stdout \
      | bcftools view \
        --apply-filters PASS,. \
        --min-ac ~{min_ac} \
        --include 'INFO/AN >= ~{min_an}' \
        --samples-file ~{samples_list} \
        --force-samples \
      | bcftools +fill-tags -- -t AC,AN,AF,F_MISSING \
      | bcftools view \
        --min-ac ~{min_ac} \
        --min-af $global_min_af \
        --max-af $global_max_af \
        --include 'INFO/AN >= ~{min_an} & INFO/F_MISSING <= ~{max_ncr}' \
      | bcftools annotate -x ^FORMAT/GT,FORMAT/AD,FORMAT/DP \
        -Oz -o local_vcfs/$idx.vcf.gz
      tabix -p vcf local_vcfs/$idx.vcf.gz

      # Delete VCF & index if no variants are found
      if [ $( bcftools query -f '%CHROM\n' local_vcfs/$idx.vcf.gz | wc -l ) -eq 0 ]; then
        # Except for the first shard, which we save in case we need a dummy header-only VCF later
        if [ $idx -eq 1 ]; then
          cp local_vcfs/$idx.vcf.gz dummy.vcf.gz
          cp local_vcfs/$idx.vcf.gz.tbi dummy.vcf.gz.tbi
        fi
        rm local_vcfs/$idx.vcf.gz local_vcfs/$idx.vcf.gz.tbi
      fi

    done < <( awk -v FS="\t" -v OFS="\t" '{ print NR, $1, $2 }' ~{snv_info_tsv} )
    find local_vcfs/ -name "*.vcf.gz" \
    | sort -V | awk -v OFS="\t" '{ print NR, $1 }' \
    > local_vcfs.list

    # Check to ensure some VCFs overlapped query intervals; otherwise, make dummy VCF and exit
    if [ $( cat local_vcfs.list | wc -l ) -eq 0 ]; then

      cp dummy.vcf.gz ~{out_vcf}
      cp dummy.vcf.gz.tbi ~{out_tbi}
    
    else

      # Combine all local mini-VCFs per query region, with more precise filtering
      mkdir region_vcfs
      while read ridx chrom start end minaf maxaf; do
        
        echo -e "\nMerging local VCFs for region $ridx ($chrom:$start-$end)..."
        
        while read vidx vcf; do

          bcftools view \
            --regions "$chrom:$start-$end" \
            --min-af $minaf \
            --max-af $maxaf \
            -Oz -o region_vcfs/$ridx.$vidx.vcf.gz \
            $vcf
          tabix -p vcf -f region_vcfs/$ridx.$vidx.vcf.gz

        done < local_vcfs.list

        find region_vcfs/ -name "$ridx.*.vcf.gz" \
        | sort -V > region_vcfs/$ridx.input_vcfs.list
        
        bcftools concat \
          --threads ~{concat_threads} \
          --allow-overlaps \
          --remove-duplicates \
          --file-list region_vcfs/$ridx.input_vcfs.list \
          -Oz -o tmp.concat.vcf.gz
        bcftools sort \
          --max-mem "~{sort_mem_mb}M" \
          --temp-dir region_vcfs/ \
          -Oz -o region_vcfs/cleaned.$ridx.vcf.gz \
          tmp.concat.vcf.gz
        tabix -p vcf region_vcfs/cleaned.$ridx.vcf.gz
        rm tmp.concat.vcf.gz

        rm region_vcfs/$ridx.*.vcf.gz* region_vcfs/$ridx.input_vcfs.list

      done < ~{query_intervals}

      # Finally, combine cleaned & sorted VCFs across all regions
      find region_vcfs/ -name "cleaned.*.vcf.gz" \
      | sort -V > outer_merge_inputs.list
      bcftools concat \
        --threads ~{concat_threads} \
        --naive \
        --file-list outer_merge_inputs.list \
        -Oz -o ~{out_vcf}
      tabix -p vcf -f ~{out_vcf}

    fi
  >>>

  output {
    File snv_vcf = out_vcf
    File snv_vcf_idx = out_tbi
  }

  runtime {
    docker: gatk_docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb +" HDD"
    preemptible: n_preemptible
    maxRetries: 1
  }
}


# Bifurcates an SV VCF based on frequency and number of alleles
task SplitSvs {
  input {
    File vcf
    File vcf_idx

    Float min_af
    Float max_af
    Int min_ac
    Int max_ac
    Int min_an

    String output_prefix

    String g2c_analysis_docker
  }

  String elig_outfile = output_prefix + ".regeno_eligible_svs.vcf.gz"
  String pt_outfile = output_prefix + ".passthrough_svs.vcf.gz"
  Int disk_gb = (3 * ceil(size([vcf], "GB"))) + 10

  command <<<
    set -eu -o pipefail

    /opt/pancan_germline_wgs/scripts/variant_filtering/bifurcate_svs_for_regenotyping.py \
      -i ~{vcf} \
      -e "~{elig_outfile}" \
      -p "~{pt_outfile}" \
      --min-af ~{min_af} \
      --max-af ~{max_af} \
      --min-ac ~{min_ac} \
      --max-ac ~{max_ac} \
      --min-an ~{min_an}

    tabix -p vcf "~{elig_outfile}"
    tabix -p vcf "~{pt_outfile}"
  >>>

  output {
    File target_sv_vcf = "~{elig_outfile}"
    File target_sv_vcf_idx = "~{elig_outfile}.tbi"
    File passthrough_sv_vcf = "~{pt_outfile}"
    File passthrough_sv_vcf_idx = "~{pt_outfile}.tbi"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "3.75 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 1
    maxRetries: 1
  }
}


# Writes a two-column .tsv of SNV VCF and index information
task WriteSnvInfo {
  input {
    Array[String] vcf_uris
    Array[String] tbi_uris
    String output_prefix
    String docker
  }

  String outfile = output_prefix + ".snv_vcf_info.tsv"

  command <<<
    set -eu -o pipefail

    paste \
      ~{write_lines(vcf_uris)} \
      ~{write_lines(tbi_uris)} \
    > ~{outfile}
  >>>

  output {
    File snv_info_tsv = outfile
  }

  runtime {
    docker: docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk 25 HDD"
    preemptible: 1
    maxRetries: 1
  }
}
