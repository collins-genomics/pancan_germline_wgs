# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Common WDL tasks shared across workflows


version 1.0


import "https://raw.githubusercontent.com/collins-genomics/pancan_germline_wgs/main/wdl/Utilities.wdl" as Tasks


workflow ConcatVcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    String out_prefix

    String bcftools_concat_options = ""

    Float mem_gb = 3.5
    Int cpu_cores = 2
    Int? disk_gb

    String bcftools_docker
  }  

  call Tasks.ConcatVcfs as Concat {
    input:
      vcfs = vcfs,
      vcf_idxs = vcf_idxs,
      out_prefix = out_prefix,
      bcftools_concat_options = bcftools_concat_options,
      mem_gb = mem_gb,
      cpu_cores = cpu_cores,
      disk_gb = disk_gb,
      bcftools_docker = bcftools_docker
  }

  output {
    File output_vcf = Concat.merged_vcf
    File output_vcf_idx = Concat.merged_vcf_idx
  }
}