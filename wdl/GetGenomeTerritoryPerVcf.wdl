# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Compute chromosomal regions spanned by one or more VCFs


version 1.0


import "Utilities.wdl" as Utils


workflow GetGenomeTerritoryPerVcf {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs

    Boolean compute_density = true  # Should territory density be computed? (Can be useful for checking for overlapping territories)
    File? genome_file               # BEDTools-style genome file. Required if compute_density is true

    String output_prefix

    String g2c_analysis_docker
  }

  # Scatter over VCFs and process in parallel
  scatter ( vcf_info in zip(vcfs, vcf_idxs) ) {
    call GetTerritory {
      input:
        vcf = vcf_info.left,
        vcf_idx = vcf_info.right,
        bcftools_docker = g2c_analysis_docker
    }
  }

  # Compute density map of intervals across all VCFs
  if ( compute_density && defined(genome_file) ) {
    call Utils.CalcBedDensity as CalcDensity {
      input:
        beds = GetTerritory.territory_bed,
        genome_file = select_first(select_all([genome_file])),
        output_prefix = output_prefix,
        bedtools_docker = g2c_analysis_docker
    }
  }
  if ( !compute_density || !defined(genome_file) ) {
    call Utils.ConcatTextFiles as ConcatTerritories {
      input:
        shards = GetTerritory.territory_bed,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = false,
        output_filename = output_prefix + ".territories.bed.gz",
        docker = g2c_analysis_docker
    }
  }
  
  output {
    File territories = select_first([CalcDensity.density_bed, ConcatTerritories.merged_file])
  }
}


task GetTerritory {
  input {
    File vcf
    File vcf_idx

    String bcftools_docker
  }

  String outfile = basename(vcf, ".vcf.gz") + ".territories.bed.gz"
  Int disk_gb = ceil(2 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    bcftools query -f '%CHROM\t%POS\n' ~{vcf} \
    | awk -v OFS="\t" \
      '!($1 in first) { first[$1] = $2 }
      { last[$1] = $2 }
      END {
        for (chrom in first)
          print chrom, first[chrom], last[chrom]
      }' | sort -Vk1,1 -k2,2n -k3,3n \
    | bedtools merge -i - \
    | bgzip -c \
    > "~{outfile}"
  >>>

  output {
    File territory_bed = "~{outfile}"
  }

  runtime {
    docker: bcftools_docker
    memory: "3.75 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 1
    maxRetries: 1
  }
}