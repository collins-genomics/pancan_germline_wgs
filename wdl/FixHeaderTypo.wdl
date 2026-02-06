# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Fix stupid typo in post-imputation SV VCF header


version 1.0


workflow FixHeaderTypo {
  input {
    File vcf
    File vcf_idx
    String bcftools_docker
  }

  call FixTypo {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      docker = bcftools_docker
  }

  output {
    File fixed_vcf = FixTypo.fixed_vcf
    File fixed_vcf_idx = FixTypo.fixed_vcf_idx
  }
}


task FixTypo {
  input {
    File vcf
    File vcf_idx
    String docker
  }

  Int disk_gb = ceil(3 * size(vcf, "GB")) + 10
  String outfile = basename(vcf, ".vcf.gz") + ".typo_fixed.vcf.gz"

  command <<<
    set -eu -o pipefail

    bcftools view --header-only ~{vcf} \
    | sed 's/Decription/Description/g' \
    > header.vcf

    bcftools reheader -h header.vcf -o ~{outfile} ~{vcf}
    tabix -p vcf -f ~{outfile}
  >>>

  output {
    File fixed_vcf = outfile
    File fixed_vcf_idx = "~{outfile}.tbi"
  }


  runtime {
    docker: docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 1
    maxRetries: 1
  }
}