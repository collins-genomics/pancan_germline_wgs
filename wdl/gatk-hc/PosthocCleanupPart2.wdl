# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Post hoc outlier sample exclusion and site-level hard filtering for GATK-HC output

# Part 1 of 2, as follows:
# 1. Variant normalization & basic site-level hard filters
# [Manual step between] Outlier identification
# 2. Outlier sample exclusion & reapplication of basic hard filters


version 1.0


import "Utilities.wdl" as Utils


workflow PosthocCleanupPart2 {
  input {
    # Input can be specified in two ways: 
    # 1. Arrays of VCFs & tabix indexes:
    Array[File]? vcfs
    Array[File]? vcf_idxs
    # 2. Or a two-column .tsv of VCF and index GCS URIs
    File? vcf_info_tsv

    File exclude_samples_list

    # Parallelization control (makes workflow metadata easier to query)
    Int vcfs_per_shard = 10

    String linux_docker
    String g2c_analysis_docker
  }

  
  # If inputs are defined as arrays, write as .tsv for subsequent chunking
  if ( ! defined(vcf_info_tsv) ) {
    call Utils.VcfArrayToTsv {
      input:
        vcfs = select_first([vcfs, []]),
        vcf_idxs = select_first([vcf_idxs, []]),
        output_prefix = "PosthocCleanupPart2",
        docker = linux_docker
    }
  }

  # Shard VCF URI list
  call Utils.ShardTextFile as ShardVcfList {
    input:
      input_file = select_first([vcf_info_tsv, VcfArrayToTsv.vcf_info_tsv]),
      lines_per_split = vcfs_per_shard,
      out_prefix = "PosthocCleanupPart2.",
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

      File in_vcf_shard = inner_vcf_infos[j].left
      File in_vcf_shard_idx = inner_vcf_infos[j].right

      call CleanupPart2 {
        input:
          vcf = in_vcf_shard,
          vcf_idx = in_vcf_shard_idx,
          exclude_samples_list = exclude_samples_list,
          docker = g2c_analysis_docker
      }
    }
  }

  output {
    Array[File] filtered_vcfs = flatten(CleanupPart2.filtered_vcf)
    Array[File] filtered_vcf_idxs = flatten(CleanupPart2.filtered_vcf_idx)
  }
}


task CleanupPart2 {
  input {
    File vcf
    File vcf_idx
    File exclude_samples_list

    String docker
    Float mem_gb = 3.75
    Int n_cpu = 2
  }

  Int disk_gb = ceil(2.25 * size([vcf], "GB")) + 10
  String outfile = basename(vcf, ".vcf.gz") + ".posthoc_filtered.vcf.gz"

  command <<<
    set -euo pipefail

    # Post hoc filters as follows:
    # 1. Exclude outlier samples
    # 2. Exclude variants with no non-reference genotypes with DP > 10 and GQ > 20

    bcftools view \
      --samples-file "^~{exclude_samples_list}" \
      ~{vcf} \
    | bcftools view \
      --include 'INFO/AC>0 & FORMAT/DP>10 & FORMAT/GQ>20 & GT="alt" & alt[0] != "*"' \
      -Oz -o "~{outfile}"

    tabix -p vcf "~{outfile}"
  >>>

  output {
    File filtered_vcf = "~{outfile}"
    File filtered_vcf_idx = "~{outfile}.tbi"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    max_retries: 1
  }
}

