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
    Int min_snv_ac = 1                   # Minimum SNV AC to be considered. Applied in addition to snv_freq_scalar criteria.
    File? snv_exclusion_bed              # SNVs overlapping this BED file will be excluded

    # Imputation parameters
    Int max_snps_per_flank = 10          # Max number of SNPs to include per flank
    Float min_ld_r2 = 0.2                # Minimum LD between SNV & SV to permit for imputation
    String ref_build = "hg38"            # Plink-styled reference indicator
    Float min_carrier_accuracy = 0.7     # Minimum (carrier|ref) accuracy to accept an SV imputation model as well-fit
    Float min_imputation_r2 = 0.3        # Minimum R2 between raw & imputed SV allele dosages to accept model as well-fit
    File? training_samples_list          # Optional list of sample IDs to consider when computing LD and training imputation models
    File? sample_group_labels            # Optional two-column .tsv mapping sample IDs to major group labels (e.g., continental ancestry). No header.
    File? sample_covariates              # Optional .tsv of sample IDs + any technical covariates for SV imputation. Has header.

    # Parallelization options
    Int svs_per_shard = 200
    Int snv_vcfs_per_shard = 100

    File genome_file                     # BEDTools-style .genome file

    String output_prefix

    String g2c_analysis_docker
    String linux_docker
  }

  # Determine method of SNV VCF input
  if ( defined(snv_vcfs) && defined(snv_vcf_idxs) ) {
    call Utils.WriteVcfInfo as WriteSnvInfo {
      input:
        vcf_uris = select_first([snv_vcfs, []]),
        tbi_uris = select_first([snv_vcf_idxs, []]),
        output_prefix = output_prefix + ".snv",
        docker = linux_docker
    }
  }
  File snv_info_tsv = select_first(select_all([WriteSnvInfo.vcf_info_tsv, snv_vcf_info_tsv]))
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
          g2c_analysis_docker = g2c_analysis_docker
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
      if ( DefineQueryIntervals.n_intervals > 0 ) {
        Int n_chunks = length(vcf_info_chunks)
        scatter ( i in range(n_chunks) ) {
          call QuerySnvs {
            input:
              snv_info_tsv = vcf_info_chunks[i],
              query_intervals = DefineQueryIntervals.query_intervals,
              samples_list = FindSharedSamples.intersection_file,
              snv_exclusion_bed = snv_exclusion_bed,
              min_ac = select_first(select_all([min_snv_ac, floor(min_sv_ac / snv_freq_scalar)])),
              min_an = min_an,
              max_ncr = 1 - min_snv_call_rate,
              output_prefix = shard_prefix + ".chunk_" + i
          }
        }
        if ( n_chunks > 1 ) {
          call Utils.ConcatVcfs as ConcatSnvs {
            input:
              vcfs = QuerySnvs.snv_vcf,
              vcf_idxs = QuerySnvs.snv_vcf_idx,
              out_prefix = shard_prefix + ".eligible_snvs",
              bcftools_concat_options = "--allow-overlaps --remove-duplicates",
              check_index_localization = false,
              mem_gb = 1.75,
              bcftools_docker = g2c_analysis_docker
          }
        }
        File merged_snv_vcf = select_first(flatten([[ConcatSnvs.merged_vcf], QuerySnvs.snv_vcf]))
        File merged_snv_vcf_idx = select_first(flatten([[ConcatSnvs.merged_vcf_idx], QuerySnvs.snv_vcf_idx]))

        # Define subset of samples to use for LD computation and imputation model training, if optioned
        if ( defined(training_samples_list) ) {
          Array[File] training_sample_lists_preint = select_all([FindSharedSamples.intersection_file, training_samples_list])
          call Utils.IntersectTextFiles as DefineTrainingSamples {
            input:
              files = training_sample_lists_preint,
              outfile = output_prefix + ".training_samples.list",
              docker = linux_docker
          }
        }

        # Compute LD for each SV, extract AD matrixes, fit regression model, and predict GTs for all samples
        call ImputeSvs {
          input:
            sv_vcf = sv_training_vcf,
            sv_vcf_idx = sv_training_vcf_idx,
            snv_vcf = merged_snv_vcf,
            snv_vcf_idx = merged_snv_vcf_idx,
            inclusion_samples_list = FindSharedSamples.intersection_file,
            training_samples_list = DefineTrainingSamples.intersection_file,
            sample_group_labels = sample_group_labels,
            sample_covariates = sample_covariates,
            breakpoint_buffer_bp = breakpoint_buffer_bp,
            breakpoint_window_bp = breakpoint_window_bp,
            snv_freq_scalar = snv_freq_scalar,
            min_snv_ac = min_snv_ac,
            min_ld_r2 = min_ld_r2,
            min_accuracy = min_carrier_accuracy,
            min_imputation_r2 = min_imputation_r2,
            min_sv_ac = min_sv_ac,
            ref_build = ref_build,
            max_snps_per_flank = max_snps_per_flank,
            mask_training_sv_gts = mask_training_sv_gts,
            sv_mask_field = sv_training_mask_field,
            output_prefix = shard_prefix,
            g2c_analysis_docker = g2c_analysis_docker
        }

        # Update SV GTs with imputed results
        if ( ImputeSvs.n_imputed_svs > 0 ) {
          call UpdateGts {
            input:
              vcf = vcf_info.left,
              vcf_idx = vcf_info.right,
              updates_tsv = ImputeSvs.imputation_results,
              g2c_analysis_docker = g2c_analysis_docker
          }
        }
      }
    }

    File imputed_vcf = select_first([UpdateGts.updated_vcf, vcf_info.left])
    File imputed_vcf_idx = select_first([UpdateGts.updated_vcf_idx, vcf_info.right])
  }

  # Concatenate all updated SV VCFs with the passthrough VCF
  call Utils.ConcatVcfs {
    input:
      vcfs = flatten([imputed_vcf, [SplitSvs.passthrough_sv_vcf]]),
      vcf_idxs = flatten([imputed_vcf_idx, [SplitSvs.passthrough_sv_vcf_idx]]),
      out_prefix = output_prefix + ".imputed",
      bcftools_concat_options = "-a -D",
      check_index_localization = true,
      bcftools_docker = g2c_analysis_docker
  }

  # Concatenate & compress all logs for archival
  Array[File] all_logs = select_all(ImputeSvs.imputation_log)
  if ( length(all_logs) > 0 ){
    call Utils.ConcatTextFiles as ConcatLogs {
      input:
        shards = all_logs,
        compression_command = "gzip -c",
        output_filename = output_prefix + ".imputation_logs.tsv.gz",
        docker = linux_docker
    }
  }

  output {
    File refined_vcf = ConcatVcfs.merged_vcf
    File refined_vcf_idx = ConcatVcfs.merged_vcf_idx
    File? imputation_logs = ConcatLogs.merged_file
  }
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

    cat ~{outfile} | wc -l > n_intervals.txt
  >>>

  output {
    File query_intervals = outfile
    Int n_intervals = read_int("n_intervals.txt")
  }

  runtime {
    docker: bcftools_docker
    memory: "1.75 GB"
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
    
    File inclusion_samples_list
    File? training_samples_list
    File? sample_group_labels
    File? sample_covariates
    Int breakpoint_buffer_bp
    Int breakpoint_window_bp
    Float snv_freq_scalar

    Float min_ld_r2
    Int min_sv_ac
    Int min_snv_ac = 1
    Float min_accuracy
    Float min_imputation_r2
    String ref_build
    Int max_snps_per_flank

    Boolean mask_training_sv_gts
    String sv_mask_field
    Int sv_mask_retries = 2
    
    String output_prefix

    String g2c_analysis_docker

    Float mem_gb = 3.5
    Int n_cpu = 2
    Int n_preemptible = 1
  }

  File training_samples_list_use = if defined(training_samples_list) then training_samples_list else inclusion_samples_list

  String groups_cmd = if defined(sample_group_labels) then "--sample-group-labels " + basename(select_first([sample_group_labels, ""])) else ""
  String covars_cmd = if defined(sample_covariates) then "--sample-covariates " + basename(select_first([sample_covariates, ""])) else ""

  Int min_snv_ac_nofloor = floor(min_sv_ac / snv_freq_scalar)
  Int min_snv_ac_use = if min_snv_ac_nofloor < min_snv_ac then min_snv_ac else min_snv_ac_nofloor

  Array[Int] retry_counter = if mask_training_sv_gts then range(sv_mask_retries + 1) else [0]

  String outfile = output_prefix + ".imputation_results.tsv.gz"
  String out_log  = output_prefix + ".sv_imputation.log"

  Int disk_gb = ceil(2 * size([sv_vcf, snv_vcf], "GB")) + 10

  Int bcftools_threads = (2 * n_cpu) - 1

  command <<<
    set -eu -o pipefail

    # Make list of SVs to process
    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\t%INFO/AF\n' ~{sv_vcf} \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
    > svs.bed

    # Make dummy .fam file
    awk -v OFS="\t" '{ print $1, $1, 0, 0, 0, 0 }' ~{inclusion_samples_list} > samples.fam

    # Localize optional files
    if ~{defined(sample_group_labels)}; then
      mv ~{sample_group_labels} ./
    fi
    if ~{defined(sample_covariates)}; then
      mv ~{sample_covariates} ./
    fi

    # Make central directory for holding all final imputation results
    mkdir imp_res
    echo -e "#sv_id\tsample\tGT\tGQ\tAD" > imp_res_header.tsv

    # If no eligible SVs are found, write empty outputs and exit early
    if [ $( cat svs.bed | wc -l ) -eq 0 ]; then
      echo "No eligible SVs found. Writing empty output and exiting."
      cat imp_res_header.tsv | gzip -c > ~{outfile}
      echo "0" > imputed_svs.count.txt
      touch ~{out_log}
      exit 0
    fi

    # Process each SV in serial
    while read chrom start end svid svaf; do

      echo -e "\nNow processing $svid..."

      # AF must be defined for variant to be processed
      if [ -z $svaf ]; then
        echo -e "SV AF undefined for $svid; skipping..."
        continue
      fi

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
      min_snv_af=$( echo $svaf | awk -v scalar=~{snv_freq_scalar} '{ af=($1 / scalar); if (af<0) af=0; print af }' )
      max_snv_af=$( echo $svaf | awk -v scalar=~{snv_freq_scalar} '{ af=($1 * scalar); if (af>1) af=1; print af }' )

      # Extract SNPs and SV record
      for flank in left right; do
        bcftools view \
          --samples-file ~{inclusion_samples_list} \
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
        --samples-file ~{inclusion_samples_list} \
        --force-samples \
        --include "ID = \"$svid\"" \
        -Oz -o $svid/sv.vcf.gz \
        ~{sv_vcf}
      tabix -p vcf -f $svid/sv.vcf.gz

      # Enter retry loop for increasingly strict SV GT filtering
      for retry_k in ~{sep=" " retry_counter}; do

        # Apply SV GT filter if not the first try
        masked_ac=$( bcftools query -f '%INFO/AC\n' $svid/sv.vcf.gz )
        if [ $retry_k -gt 0 ]; then
          echo -e "\nRetry number $retry_k of stricter SV GT masking..."
          /opt/pancan_germline_wgs/scripts/variant_filtering/mask_sv_gts_for_regenotyping.py \
            --input-vcf $svid/sv.vcf.gz \
            --quality-field "~{sv_mask_field}" \
          | bcftools annotate -x INFO/AC,INFO/AN,INFO/AC \
          | bcftools +fill-tags \
            -Oz -o $svid/sv.masked.vcf.gz \
            -- -t AC,AN,AF
          masked_ac=$( bcftools query -f '%INFO/AC\n' $svid/sv.masked.vcf.gz )
          echo -e "Stricter SV GT masking retained $masked_ac non-reference alleles..."
          if [ $masked_ac -ge ~{min_sv_ac} ]; then
            mv $svid/sv.masked.vcf.gz $svid/sv.vcf.gz
            tabix -p vcf -f $svid/sv.vcf.gz
          fi
        fi
        if [ $masked_ac -lt ~{min_sv_ac} ]; then
          touch $svid/ld.vcor.slim
          break
        fi

        # Make VCF sandwich of (left flanking SNVs) + SV + (right flanking SNVs)
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
          --ld-window-kb $(( 2 * ( ~{breakpoint_buffer_bp} + ~{breakpoint_window_bp} + $svlen + 1 ) / 1000 )) \
          --ld-window-r2 ~{min_ld_r2} \
          --split-par "~{ref_build}" \
          --polyploid-mode missing \
          --fam samples.fam \
          --ld-snp "$svid" \
          --vcf $svid/sandwich.vcf.gz \
          --out $svid/ld
        cut -f6,7 $svid/ld.vcor | sed '1d' | sort -nrk2,2 -k1,1V > $svid/ld.vcor.slim

        # If no tag SNPs are found, do nothing more
        if [ $( cat $svid/ld.vcor.slim | wc -l ) -eq 0 ]; then
          cp imp_res_header.tsv imp_res/$svid.imputation_results.tsv
          echo -e "Found no tag SNPs for $svid; continuing...\n"

        # Otherwise, continue with imputation
        else

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

          # Impute SV GTs
          echo -e "Now imputing $svid..."
          /opt/pancan_germline_wgs/scripts/variant_filtering/impute_sv_gts.R \
            --ad $svid/$svid.ad.tsv.gz \
            --sv-id "$svid" \
            ~{covars_cmd} \
            ~{groups_cmd} \
            --training-samples ~{training_samples_list_use} \
            --min-ac ~{min_sv_ac} \
            --min-snv-ac ~{min_snv_ac_use} \
            --min-accuracy ~{min_accuracy} \
            --min-r2 ~{min_imputation_r2} \
            --out-tsv imp_res/$svid.imputation_results.tsv \
          >> ~{out_log} || true

        fi

        if [ -s imp_res/$svid.imputation_results.tsv ] && \
           [ $( cat imp_res/$svid.imputation_results.tsv | wc -l ) -gt 1 ]; then
          break
        fi

      done

      # Clean up
      rm -rf $svid

    done < svs.bed

    # Concatenate imputation results
    cat imp_res/*.tsv \
    | { grep -ve '^#' || true; } \
    | cat imp_res_header.tsv - \
    | gzip -c > ~{outfile}
    n_imputed=$( zcat ~{outfile} | grep -ve '^#' | cut -f1 \
                 | sort | uniq | wc -l || true )
    echo "$n_imputed" > imputed_svs.count.txt
  >>>

  output {
    File imputation_results = outfile
    Int n_imputed_svs = read_int("imputed_svs.count.txt")
    File? imputation_log = out_log
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

    Float mem_gb = 1.75

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
    memory: "~{mem_gb} GB"
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
    File? samples_list
    
    Int min_ac
    Int min_an
    Float max_ncr
    File? snv_exclusion_bed

    Int max_stream_attempts = 5

    String output_prefix

    Int disk_gb = 30
    Float mem_gb = 7.5
    Int n_cpu = 2
    Int n_preemptible = 1
    String gatk_docker = "broadinstitute/gatk:4.6.2.0"
  }

  String out_vcf = output_prefix + ".snvs.vcf.gz"
  String out_tbi = out_vcf + ".tbi"

  String samples_cmd = if defined(samples_list) 
                       then "--force-samples --samples-file " + basename(select_first([samples_list, ""]))
                       else ""

  String excl_cmd = if defined(snv_exclusion_bed) 
                    then "--targets-file ^" + basename(select_first([snv_exclusion_bed, ""]))
                    else ""

  Int concat_threads = 2 * n_cpu
  Int sort_mem_mb = floor(1000 * (mem_gb / 3))

  command <<<
    set -eu -o pipefail

    # Relocate sample list and exclusion bed, if provided
    if ~{defined(samples_list)}; then
      cp ~{default=" " samples_list} ./
    fi
    if ~{defined(snv_exclusion_bed)}; then
      cp ~{default=" " snv_exclusion_bed} ./
    fi

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
        ~{samples_cmd} \
        ~{excl_cmd} \
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


# Injects imputed GTs into an SV VCF
task UpdateGts {
  input {
    File vcf
    File vcf_idx
    File updates_tsv

    Int gq_offset = 9

    String g2c_analysis_docker
  }

  String outfile = basename(vcf, ".vcf.gz") + ".imputed.vcf.gz"

  Int disk_gb = ceil(3 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    # Gather original GT information for all samples for qualifying records
    zcat ~{updates_tsv} | cut -f1 | grep -ve '^#' | sort -V | uniq > svids.list
    echo -e "#sv_id\tsample\tGT\tGQ" > current.header
    bcftools query \
      -i 'ID=@svids.list' \
      -f '[%ID\t%SAMPLE\t%GT\t%GQ\n]' \
      ~{vcf} \
    | sort -Vk1,1 -k2,2V \
    | cat current.header - \
    | gzip -c \
    > current_gts.tsv.gz

    # Array-based comparison of old & imputed GTs to prioritize GTs to update
    # (Vastly faster than doing this serially in pysam)
    /opt/pancan_germline_wgs/scripts/variant_filtering/filter_imputed_gts.R \
      --old-gts current_gts.tsv.gz \
      --imputed-gts ~{updates_tsv} \
      --gq-offset ~{gq_offset} \
      --out-tsv filtered.updates.tsv
    gzip -f filtered.updates.tsv

    # Only update the GTs passing the above filter
    /opt/pancan_germline_wgs/scripts/variant_filtering/inject_imputed_sv_gts.py \
      -i ~{vcf} \
      -u filtered.updates.tsv.gz \
      --no-compare \
    | bcftools +fill-tags -- -t AC,AN,AF \
    | bcftools view \
      -i 'AC > 0 | FILTER = "MULTIALLELIC"' \
      -Oz -o ~{outfile}
    tabix -p vcf -f ~{outfile}
  >>>

  output {
    File updated_vcf = outfile
    File updated_vcf_idx = outfile + ".tbi"
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

