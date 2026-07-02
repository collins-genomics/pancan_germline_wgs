# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Collect features for genotype filtering model for a single input VCF
# Given the design of G2C, this workflow assumes input VCF contains 
# exactly one variant class (SNV, indel, SV)


version 1.0


workflow CollectGTFilterFeatures {
	input {
    File vcf
    File vcf_idx

    String g2c_analysis_docker
  }

  # Collect no-call rates per variant class
  call CollectNoCallRates {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Collect site-level features
  # TODO: implement this

  # Collect genotype-level features
  # TODO: implement this

  output {
    File nocall_counts = CollectNoCallRates.nocall_counts
  }
}


# Compute no-call rates for all samples
task CollectNoCallRates {
  input {
    File vcf
    File vcf_idx

    String g2c_analysis_docker
  }

  String outfile = basename(vcf, ".vcf.gz") + ".nocall_counts.tsv.gz"
  Int disk_gb = ceil(1.5 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    # Get list of all samples present in VCF
    bcftools query -l ~{vcf} > samples.list

    # Get total number of records in file
    n_records=$( bcftools index -n ~{vcf} )

    # Tabulate missing genotypes
    bcftools query \
      -i 'INFO/SVTYPE != "MCNV" & GT="mis"' \
      -f '[%SAMPLE\n]' \
      ~{vcf} \
    | cat samples.list - \
    | sort -Vk1,1 \
    | uniq -c \
    | awk -v n="$n_records" -v OFS="\t" '{ print $2, n, $1-1 }' \
    | gzip -c \
    > "~{outfile}"
  >>>

  output {
    File nocall_counts = outfile
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    maxRetries: 1
  }
}
