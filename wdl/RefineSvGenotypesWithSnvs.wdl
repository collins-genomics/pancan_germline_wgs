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

    # SV/SNV filtering parameters
    Float min_sv_af = 0.05
    Int min_sv_ac = 20
    Int breakpoint_buffer_bp = 5000
    Int breakpoint_window_bp = 100000
    File? snv_exclusion_bed

    Int svs_per_shard = 100

    String output_prefix

    String g2c_pipeline_docker
    String linux_docker
  }

  # Determine method of SNV VCF input
  if ( defined(snv_vcfs) && defined(snv_vcf_idxs) ) {
    call WriteSnvInfo {
      input:
        vcf_uris = snv_vcfs,
        tbi_uris = snv_vcf_idxs,
        output_prefix = output_prefix,
        docker = linux_docker
    }
  }
  File snv_info_tsv = select_first(select_all([WriteSnvInfo.snv_info_tsv, snv_vcf_info_tsv]))
  if ( !defined(snv_vcfs) ) {
    call GetFirstSnvVcf {
      input:
        tsv = snv_vcf_info_tsv,
        docker = linux_docker
    }
  }
  File first_snv_vcf = select_first(flatten(select_all([select_all(snv_vcfs), 
                                            select_all([GetFirstSnvVcf.vcf_uri])])))
  File first_snv_idx = select_first(flatten(select_all([select_all(snv_vcf_idxs), 
                                            select_all([GetFirstSnvVcf.tbi_uri])])))

  # Gather list of matching samples between snv_vcfs and sv_vcf
  call Utils.StreamSamplesFromVcfHeader as GetSvSamples {
    input:
      vcf = sv_vcf,
      vcf_idx = sv_vcf_idx,
      bcftools_docker = g2c_pipeline_docker
  }
  call Utils.StreamSamplesFromVcfHeader as GetSnvSamples {
    input:
      vcf = first_snv_vcf,
      vcf_idx = first_snv_idx,
      bcftools_docker = g2c_pipeline_docker
  }
  call Utils.IntersectTextFiles as FindSharedSamples {
    input:
      files = [GetSvSamples.sample_list, GetSnvSamples.sample_list],
      outfile = output_prefix + ".shared_samples.list",
      docker = linux_docker
  }

  # Separate qualifying and non-qualifying SVs
  call SplitSvs {
    input:
      vcf = sv_vcf,
      vcf_idx = sv_vcf_idx,
      min_af = min_sv_af,
      max_af = 1 - min_sv_af,
      min_ac = min_sv_ac,
      max_ac = (2 * GetSvSamples.n_samples) - min_sv_ac,
      output_prefix = basename(sv_vcf, ".vcf.gz"),
      g2c_pipeline_docker = g2c_pipeline_docker
  }

  # Shard qualifying SV VCF for parallel processing
  call Utils.ShardVcf as ShardTargetSvs {
    input:
      vcf = SplitSvs.target_sv_vcf,
      vcf_idx = SplitSvs.target_sv_vcf_idx,
      records_per_shard = svs_per_shard,
      bcftools_docker = g2c_pipeline_docker,
      n_preemptible = 1
  }

  # Process each qualifying SV VCF in parallel
  scatter ( vcf_info in zip(ShardVcf.vcf_shards, ShardVcf.vcf_shard_idxs) ) {

    # Filter SNVs
    call QuerySnvs {
      input:
        sv_vcf = vcf_info.left,
        sv_vcf_idx = vcf_info.right,
        snv_info_tsv = snv_info_tsv,
        breakpoint_buffer_bp = breakpoint_buffer_bp,
        breakpoint_window_bp = breakpoint_window_bp,
        snv_exclusion_bed = snv_exclusion_bed

    }
    # TODO: implement this
    # - Query SNVs to the left of POS and right of END (buffer in windows of 5kb-100kb away)
    #   - Filter on min call rate
    #   - Mask low-complexity regions & segdups
    #   - Filter on AF ~ [1/5, 5] * SV_AF
    #   - Biallelic filter PASS
    # Should be clever about this -- can maybe make bedgraph or bed with min AF / min AC for any relevant SV for that window?

    # Compute LD for each SV, extract AD matrixes, fit regression model, and predict GTs for all samples
    # TODO: implement this
    # - Make VCF sandwich of (left flanking SNVs) + SV + (right flanking SNVs)
    # - Compute all LD with plink (min R2 > 0.2?)
    # - Rank-order SNVs by LD R2 per flank and take up to 5 best tag SNVs from each flank
    # - Extract allele dosage for each SNP (2 * AB)
    # - Load tag SNP AD and SV GTs into R
    # - Fit linear regression of SV AC ~ tag SNP ADs using 10-fold CV
    #   - Need to think about how to handle/prespecify train/test split (by cohort etc)
    # - Predict SV AC from tag SNPs using best-fit regression model
    # - Compute SNV-based GQ by ratio of linear distances between integer AC states (or maybe multivariate gaussian)

    # Update SV GTs
    # TODO: implement this
    # - If SNV-predicted GQ > GATK-SV GQ, return (sample, SV ID, GT, GQ) to be updated in SV VCF
  }

  # Concatenate all updated SV VCFs with the passthrough VCF
  # TODO: implement this

  output {}
}


# Extracts the first SNV VCF URI from a .tsv of VCF info
task GetFirstSnvVcf {
  input {
    File tsv
    String docker
  }

  command <<<
    set -eu -o pipefail

    cat ~{read_tsv(tsv)} | sed -n '1p' > first_info.tsv
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


# Bifurcates an SV VCF based on frequency and number of alleles
task SplitSvs {
  input {
    File vcf
    File vcf_idx

    Float min_af
    Float max_af
    Int min_ac
    Float max_ac

    String output_prefix

    String g2c_pipeline_docker
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
      --max-ac ~{max_ac}

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
    docker: g2c_pipeline_docker
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
