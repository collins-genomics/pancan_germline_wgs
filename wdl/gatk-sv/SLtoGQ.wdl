# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Updates GQ values to reflect the appropriate SL transformation
# See page 36 here:
# https://support.researchallofus.org/hc/en-us/articles/4617899955092-All-of-Us-Genomic-Quality-Report-ARCHIVED-C2022Q4R9-CDR-v7


version 1.0


import "Utilities.wdl" as Utils


workflow SLtoGQ {
  input {
    File vcf
    File vcf_idx

    Int records_per_shard = 250

    String g2c_analysis_docker
  }

  call Utils.ShardVcf {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      records_per_shard = records_per_shard,
      bcftools_docker = g2c_analysis_docker
  }

  Array[Pair[File, File]] shard_infos = zip(ShardVcf.vcf_shards, ShardVcf.vcf_shard_idxs)
  scatter ( si in shard_infos ) {
    call UpdateGq {
      input:
        vcf = si.left,
        vcf_idx = si.right,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }

  call Utils.ConcatVcfs {
    input:
      vcfs = UpdateGq.updated_vcf,
      vcf_idxs = UpdateGq.updated_vcf_idx,
      out_prefix = basename(vcf, "vcf.gz") + "gq_updated",
      bcftools_concat_options = "-a",
      bcftools_docker = g2c_analysis_docker
  }

  output {
    File updated_vcf = ConcatVcfs.merged_vcf
    File updated_vcf_idx = ConcatVcfs.merged_vcf_idx
  }
}


task UpdateGq {
  input {
    File vcf
    File vcf_idx
    String g2c_analysis_docker
  }

  String outfile = basename(vcf, "vcf.gz") + "gq_updated.vcf.gz"
  Int disk_gb = ceil(2.5 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    /opt/pancan_germline_wgs/scripts/gatksv_helpers/update_gq_from_sl.py -i "~{vcf}" -o "~{outfile}"

    tabix -p vcf -f "~{outfile}"
  >>>

  output {
    File updated_vcf = "~{outfile}"
    File updated_vcf_idx = "~{outfile}.tbi"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}
