# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Count SVs per sample per SV type
# Supports chromosome-scattered inputs for parallelization


version 1.0


import "Utilities.wdl" as Utils


workflow CountSvsPerSample {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    Int records_per_shard = 5000
    String? bcftools_view_options
    String output_prefix
    String sv_pipeline_docker
    String g2c_pipeline_docker
  }

  Array[Pair[File, File]] vcf_infos = zip(vcfs, vcf_idxs)

  scatter ( vcf_info in vcf_infos ) {
    
    if ( defined(bcftools_view_options) ) {
      call FilterVcf {
        input:
          vcf = vcf_info.left,
          vcf_idx = vcf_info.right,
          bcftools_view_options = select_first([bcftools_view_options, ""]),
          docker = sv_pipeline_docker
      }
    }

    File vcf = select_first([FilterVcf.filtered_vcf, vcf_info.left])
    File vcf_idx = select_first([FilterVcf.filtered_vcf_idx, vcf_info.right])

    call Utils.ShardVcf as ShardVcf {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx, 
        records_per_shard = records_per_shard,
        bcftools_docker = sv_pipeline_docker
    }

    Array[Pair[File, File]] shard_infos = zip(ShardVcf.vcf_shards, ShardVcf.vcf_shard_idxs)

    scatter ( shard_info in shard_infos ) {

      File shard_vcf = shard_info.left
      File shard_vcf_idx = shard_info.right

      call CountSvs {
        input:
          vcf = shard_vcf,
          vcf_idx = shard_vcf_idx,
          docker = sv_pipeline_docker
      }
    }

    if ( length(shard_infos) > 1 ) {
      call Utils.SumSvCountsPerSample as SumPerVcf {
        input:
          count_tsvs = CountSvs.counts_tsv,
          output_prefix = basename(vcf, ".vcf.gz") + ".counts.tsv",
          docker = g2c_pipeline_docker
      }
    }

    File scattered_summed_count = select_first([SumPerVcf.summed_tsv, CountSvs.counts_tsv[0]])
  }

  if ( length(vcfs) > 1 ) {
    call Utils.SumSvCountsPerSample as SumAll {
      input:
        count_tsvs = scattered_summed_count,
        output_prefix = output_prefix,
        docker = g2c_pipeline_docker
    }
  }

  output {
    File sv_counts = select_first([SumAll.summed_tsv, scattered_summed_count[0]])
  }
}


task CountSvs {
  input {
    File vcf
    File vcf_idx

    String docker
    Float mem_gb = 3.75
    Int n_cpu = 2
  }

  Int disk_gb = ceil(2 * size([vcf], "GB")) + 10
  String outfile = basename(vcf, ".vcf.gz") + ".counts.tsv"

  command <<<
    set -euo pipefail

    svtk count-svtypes ~{vcf} ~{outfile}
  >>>

  output {
    File counts_tsv = "~{outfile}"
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


task FilterVcf {
  input {
    File vcf
    File vcf_idx
    
    String bcftools_view_options

    String docker
    Float mem_gb = 3.75
    Int n_cpu = 2
  }

  Int disk_gb = ceil(2 * size([vcf], "GB")) + 10
  String outfile = basename(vcf, ".vcf.gz") + ".filtered.vcf.gz"

  command <<<
    set -euo pipefail

    bcftools view ~{bcftools_view_options} -Oz -o ~{outfile} ~{vcf}
    tabix -p vcf ~{outfile}
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
