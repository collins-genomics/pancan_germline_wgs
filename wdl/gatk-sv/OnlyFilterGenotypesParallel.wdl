# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# This is a parallel implementation of OnlyFilterGenotypes.wdl for multi-VCF input
# It skips sharding the input VCF


version 1.0


import "OnlyFilterGenotypes.wdl" as OFG


workflow OnlyFilterGenotypesParallel {
  input {
    Array[File] vcfs
    File ploidy_table

    Float no_call_rate_cutoff = 0.05  # Set to 1 to disable NCR filtering
    String sl_filter_args  # Explicitly set SL cutoffs. See apply_sl_filter.py for arguments.

    # For SanitizeHeader - bcftools annotate -x ~{header_drop_fields}
    String header_drop_fields = "FILTER/LOW_QUALITY,FORMAT/TRUTH_CN_EQUAL,FORMAT/GT_FILTER,FORMAT/CONC_ST,INFO/STATUS,INFO/TRUTH_AC,INFO/TRUTH_AN,INFO/TRUTH_AF,INFO/TRUTH_VID,INFO/CNV_CONCORDANCE,INFO/GENOTYPE_CONCORDANCE,INFO/HET_PPV,INFO/HET_SENSITIVITY,INFO/HOMVAR_PPV,INFO/HOMVAR_SENSITIVITY,INFO/MINSL,INFO/NON_REF_GENOTYPE_CONCORDANCE,INFO/SL_MAX,INFO/SL_MEAN,INFO/VAR_PPV,INFO/VAR_SENSITIVITY,INFO/VAR_SPECIFICITY"

    File? filter_script

    String sv_base_mini_docker
    String sv_pipeline_docker
  }

  scatter ( vcf in vcfs ) {
    call OFG.FilterGenotypes as FilterVcfs {
      input:
        vcf = vcf,
        output_prefix = basename(vcf, ".vcf.gz"),
        ploidy_table = ploidy_table,
        no_call_rate_cutoff = no_call_rate_cutoff,
        sl_filter_args = sl_filter_args,
        header_drop_fields = header_drop_fields,
        filter_script = filter_script,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker
    }
  }

  output {
    Array[File] filtered_vcfs = FilterVcfs.filtered_vcf
    Array[File] filtered_vcf_idxs = FilterVcfs.filtered_vcf_index
  }
}

