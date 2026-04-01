# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Wrapper to execute ApplyScriptParallelPerChrom.wdl over an array of input VCFs
# with optional concatenation


version 1.0


import "https://raw.githubusercontent.com/collins-genomics/pancan_germline_wgs/main/wdl/ApplyScriptParallelPerChrom.wdl" as Parallelize
import "https://raw.githubusercontent.com/collins-genomics/pancan_germline_wgs/main/wdl/Utilities.wdl" as Utilities


workflow MultiVcfParallelizeScript {
  input {
    Array[File] vcfs
    Array[File]? vcf_idxs

    File script
    String exec_prefix = "source"
    String script_options = ""
    Array[String?] script_files

    String? out_vcf_prefix
    Boolean scatter_by_chrom = true
    File? ref_fai
    String? inner_merge_bcftools_concat_options 

    Float? script_mem_gb
    Int? script_cpu_cores
    Int? script_disk_gb
    Float? inner_merge_mem_gb
    Int? inner_merge_cpu_cores
    Int? inner_merge_disk_gb

    Boolean merge_all_outputs = false
    String? outer_merge_bcftools_concat_options 
    Float? outer_merge_mem_gb
    Int? outer_merge_cpu_cores
    Int? outer_merge_disk_gb
    String? outer_merge_output_prefix

    String bcftools_docker
  }

  String outer_merge_prefix = select_first([outer_merge_output_prefix, basename(vcfs[0], ".vcf.gz") + "merged"])

  # Index input VCFs if necessary
  if (!defined(vcf_idxs)) {
    scatter ( vcf in vcfs ) {
      call Utilities.IndexVcf {
        input:
          vcf = vcf,
          docker = bcftools_docker
      }
    }
  }

  Array[File] definite_vcf_idxs = select_first([vcf_idxs, IndexVcf.vcf_idx])

  scatter ( vcf_info in zip(vcfs, definite_vcf_idxs) ) {

    File vcf = vcf_info.left
    File vcf_idx = vcf_info.right

    if ( scatter_by_chrom ) {
      call Parallelize.ApplyScriptParallelPerChrom as ApplyPerChrom {
        input:
          vcf = vcf,
          vcf_idx = vcf_idx,
          script = script,
          out_vcf_prefix = out_vcf_prefix,
          exec_prefix = exec_prefix,
          script_options = script_options,
          script_files = script_files,
          ref_fai = ref_fai,
          bcftools_concat_options = inner_merge_bcftools_concat_options,
          script_mem_gb = script_mem_gb,
          script_cpu_cores = script_cpu_cores,
          script_disk_gb = script_disk_gb,
          merge_mem_gb = inner_merge_mem_gb,
          merge_cpu_cores = inner_merge_cpu_cores,
          merge_disk_gb = inner_merge_disk_gb,
          bcftools_docker = bcftools_docker
      }
    }

    if ( !scatter_by_chrom ) {
      call Parallelize.ApplyScript as ApplyAsIs {
        input:
          vcf = vcf,
          vcf_idx = vcf_idx,
          script = script,
          exec_prefix = exec_prefix,
          script_options = script_options,
          script_files = script_files,
          mem_gb = script_mem_gb,
          cpu_cores = script_cpu_cores,
          disk_gb = script_disk_gb,
          docker = bcftools_docker
      }
    }
  }

  Array[File] processed_vcfs = if scatter_by_chrom then select_all(ApplyPerChrom.output_vcf) else select_all(ApplyAsIs.vcf_out)
  Array[File] processed_vcf_idxs = if scatter_by_chrom then select_all(ApplyPerChrom.output_vcf_idx) else select_all(ApplyAsIs.vcf_out_idx)

  if ( merge_all_outputs ) {
    call Utilities.ConcatVcfs {
      input:
        vcfs = processed_vcfs,
        vcf_idxs = processed_vcf_idxs,
        out_prefix = outer_merge_prefix,
        bcftools_concat_options = outer_merge_bcftools_concat_options,
        mem_gb = outer_merge_mem_gb,
        cpu_cores = outer_merge_cpu_cores,
        disk_gb = outer_merge_disk_gb,
        bcftools_docker = bcftools_docker
    }
  }

  output {
    Array[File] output_vcfs = processed_vcfs
    Array[File] output_vcf_idxs = processed_vcf_idxs
    File? merged_output_vcf = ConcatVcfs.merged_vcf
    File? merged_output_vcf_idx = ConcatVcfs.merged_vcf_idx
  }
}

