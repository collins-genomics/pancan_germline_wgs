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
    Array[File] snv_vcfs # VCFs of SNVs (can also include indels, which will be dropped)
    Array[File] snv_vcf_idxs

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

  # Gather list of matching samples between snv_vcfs and sv_vcf
  call Utils.StreamSamplesFromVcfHeader as GetSvSamples {
    input:
      vcf = sv_vcf,
      vcf_idx = sv_vcf_idx,
      bcftools_docker = g2c_pipeline_docker
  }
  call Utils.StreamSamplesFromVcfHeader as GetSnvSamples {
    input:
      vcf = snv_vcfs[0],
      vcf_idx = snv_vcf_idxs[0],
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
        snv_vcfs = snv_vcfs,
        snv_vcf_idxs = snv_vcf_idxs,
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

    # Relocate samples list to pwd, if provided
    if ~{defined(samples_list)}; then
      cp ~{default="" samples_list} ./
    fi

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
    docker: bcftools_docker
    memory: "3.75 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 1
    maxRetries: 1
  }
}


