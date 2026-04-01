# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Replacement for GATK-SV FilterBatchSamples.wdl (module 09) that excludes based on a list of sample IDs supplied by the user


version 1.0


import "https://raw.githubusercontent.com/collins-genomics/pancan_germline_wgs/main/wdl/ExcludeSamplesFromVcf.wdl" as ExcludeTask
import "https://raw.githubusercontent.com/broadinstitute/gatk-sv/v0.26.10-beta/wdl/FilterOutlierSamples.wdl" as filter_outliers
import "https://raw.githubusercontent.com/broadinstitute/gatk-sv/v0.26.10-beta/wdl/TasksMakeCohortVcf.wdl" as tasks_mcv
import "https://raw.githubusercontent.com/broadinstitute/gatk-sv/v0.26.10-beta/wdl/Utils.wdl" as util


workflow SimpleFilterBatchSamples {
   input {
    String batch
    File? manta_vcf
    File? wham_vcf
    File? melt_vcf
    File? scramble_vcf
    File? depth_vcf
    File exclude_samples_list
    String? bcftools_filter_options
    String sv_base_mini_docker
  }

  Array[File?] vcfs = [manta_vcf, wham_vcf, melt_vcf, scramble_vcf, depth_vcf]
  Array[String] algorithms = ["manta", "wham", "melt", "scramble", "depth"]  # fixed algorithms to enable File outputs to be determined
  Int num_algorithms = length(algorithms)

  scatter (i in range(num_algorithms)) {
    if (defined(vcfs[i])) {
      call ExcludeTask.ExcludeSamplesFromVcf as Exclude {
        input:
          vcf = select_first([vcfs[i]]),
          exclude_samples_list = exclude_samples_list,
          bcftools_filter_options = bcftools_filter_options,
          outfile_prefix = "~{batch}.~{algorithms[i]}.outliers_removed",
          docker = sv_base_mini_docker
      }
    }
  }

  call util.GetSampleIdsFromVcf {
    input:
      vcf = select_first(vcfs),
      sv_base_mini_docker = sv_base_mini_docker
  }

  # Write new list of samples without outliers
  call filter_outliers.FilterSampleList {
    input:
      original_samples = GetSampleIdsFromVcf.out_array,
      outlier_samples = read_lines(exclude_samples_list),
      batch = batch,
      linux_docker = sv_base_mini_docker
  }

  Array[File] pesr_vcfs_no_outliers = select_all([Exclude.filtered_vcf[0], Exclude.filtered_vcf[1], Exclude.filtered_vcf[2], Exclude.filtered_vcf[3]])
  Array[File] pesr_vcfs_no_outliers_index = select_all([Exclude.filtered_vcf_idx[0], Exclude.filtered_vcf_idx[1], Exclude.filtered_vcf_idx[2], Exclude.filtered_vcf_idx[3]])
  call tasks_mcv.ConcatVcfs as MergePesrVcfs {
    input:
      vcfs=pesr_vcfs_no_outliers,
      vcfs_idx=pesr_vcfs_no_outliers_index,
      allow_overlaps=true,
      outfile_prefix="~{batch}.filtered_pesr_merged",
      sv_base_mini_docker=sv_base_mini_docker
  }

  output {
    File? outlier_filtered_manta_vcf = Exclude.filtered_vcf[0]
    File? outlier_filtered_wham_vcf = Exclude.filtered_vcf[1]
    File? outlier_filtered_melt_vcf = Exclude.filtered_vcf[2]
    File? outlier_filtered_scramble_vcf = Exclude.filtered_vcf[3]
    File? outlier_filtered_depth_vcf = Exclude.filtered_vcf[4]

    File? outlier_filtered_manta_vcf_index = Exclude.filtered_vcf_idx[0]
    File? outlier_filtered_wham_vcf_index = Exclude.filtered_vcf_idx[1]
    File? outlier_filtered_melt_vcf_index = Exclude.filtered_vcf_idx[2]
    File? outlier_filtered_scramble_vcf_index = Exclude.filtered_vcf_idx[3]
    File? outlier_filtered_depth_vcf_index = Exclude.filtered_vcf_idx[4]

    File outlier_filtered_pesr_vcf = MergePesrVcfs.concat_vcf
    File outlier_filtered_pesr_vcf_index = MergePesrVcfs.concat_vcf_idx
    Array[String] filtered_batch_samples_list = FilterSampleList.filtered_samples_list
    File filtered_batch_samples_file = FilterSampleList.filtered_samples_file
  }
}

