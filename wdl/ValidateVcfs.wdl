# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Structurally validate one or more VCFs
# Does not check VCF spec compliance, only that the file is complete
# Corrupt or incomplete VCFs will be repaired at the last complete line


version 1.0


import "Utilities.wdl" as Utils


workflow ValidateVcfs {
  input {
    # Input can be specified in two ways: 
    # 1. Arrays of VCFs & tabix indexes:
    Array[File]? vcfs
    Array[File]? vcf_idxs
    # 2. Or a two-column .tsv of VCF and index GCS URIs
    File? vcf_info_tsv

    # Parallelization control (makes workflow metadata easier to query)
    Int vcfs_per_shard = 10

    String output_prefix = "vcf_validation"

    String linux_docker
    String g2c_analysis_docker
  }


  # If inputs are defined as arrays, write as .tsv for subsequent chunking
  if ( ! defined(vcf_info_tsv) ) {
    call Utils.VcfArrayToTsv {
      input:
        vcfs = select_first([vcfs, []]),
        vcf_idxs = select_first([vcf_idxs, []]),
        output_prefix = output_prefix,
        docker = linux_docker
    }
  }

  # Shard VCF URI list
  call Utils.ShardTextFile as ShardVcfList {
    input:
      input_file = select_first([vcf_info_tsv, VcfArrayToTsv.vcf_info_tsv]),
      lines_per_split = vcfs_per_shard,
      out_prefix = output_prefix + ".",
      g2c_analysis_docker = g2c_analysis_docker
  }
  Int n_chunks = length(ShardVcfList.shards)

  # Scatter over sharded VCF lists
  scatter ( i in range(n_chunks) ) {
    # Extract URIs from sharded file as array of strings and indexes
    call Utils.ExtractVcfArrays {
      input:
        vcf_info = ShardVcfList.shards[i],
        linux_docker = g2c_analysis_docker
    }

    Array[Pair[File, File]] inner_vcf_infos = zip(ExtractVcfArrays.vcf_uris, 
                                                  ExtractVcfArrays.vcf_tbi_uris)
    scatter ( j in range(length(inner_vcf_infos)) ) {

      File input_vcf = inner_vcf_infos[j].left
      File input_vcf_idx = inner_vcf_infos[j].right

      # Validate each VCF
      call ValidateVcf {
        input:
          vcf = input_vcf,
          vcf_idx = input_vcf_idx,
          docker = g2c_analysis_docker
      }
    }

    Array[File] valid_vcfs_inner = select_all(ValidateVcf.valid_vcf)
    Array[File] valid_vcf_idxs_inner = select_all(ValidateVcf.valid_vcf_idx)
    Array[File] repaired_vcfs_inner = select_all(ValidateVcf.repaired_vcf)
    Array[File] repaired_vcf_idxs_inner = select_all(ValidateVcf.repaired_vcf_idx)
  }

  output {
    Array[File] valid_input_vcfs = select_all(flatten(valid_vcfs_inner))
    Array[File] valid_input_vcf_idxs = select_all(flatten(valid_vcf_idxs_inner))
    Array[File] repaired_vcfs = select_all(flatten(repaired_vcfs_inner))
    Array[File] repaired_vcf_idxs = select_all(flatten(repaired_vcf_idxs_inner))
  }
}


# Validate a VCF by indexing
# If indexing fails, repairs & reindexes the VCF
# Output depends on whether the input VCF was valid
task ValidateVcf {
  input {
    File vcf
    File vcf_idx
    String docker
  }

  String vcf_basename = basename(vcf)

  String valid_vcf_path = "valid/" + vcf_basename
  String valid_vcf_idx_path = "valid/" + vcf_basename + ".tbi"
  
  String repaired_vcf_path = "repaired/" + vcf_basename
  String repaired_vcf_idx_path = "repaired/" + vcf_basename + ".tbi"
  
  Int disk_gb = ceil(3 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail
    {
      tabix -p vcf ~{vcf}
      mkdir valid
      mv ~{vcf} valid/
      mv ~{vcf_idx} valid/
    } || {
      mkdir repaired/
      /opt/pancan_germline_wgs/scripts/utilities/fix_truncated_vcf.py -i ~{vcf} \
      | bgzip -c \
      > ~{repaired_vcf_path}
      tabix -p vcf ~{repaired_vcf_path}
    }
  >>>

  output {
    File? valid_vcf = valid_vcf_path
    File? valid_vcf_idx = valid_vcf_idx_path
    File? repaired_vcf = repaired_vcf_path
    File? repaired_vcf_idx = repaired_vcf_idx_path
  }

  runtime {
    docker: docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    maxRetries: 1
  }
}

