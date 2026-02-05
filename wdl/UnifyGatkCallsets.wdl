# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Unify GATK-HC and GATK-SV callsets while collapsing large indels & small SVs


version 1.0


workflow UnifyGatkCallsets {
  input {
    Array[File] gatkhc_vcfs
    Array[File] gatkhc_vcf_idxs

    Array[File] gatksv_vcfs
    Array[File] gatksv_vcf_idxs

    File genome_file                 # BEDTools-style genome file
    File nmask_bed                   # BED file of N-masked reference intervals

    Int min_sv_size = 50
    Float size_scalar = 3            # Maximum fold-difference between sizes of indels and SVs to tolerate

    # Optional BED files to specify output sharding intervals
    # Fourth column _must_ be desired output VCF prefix
    # If not provided, outputs will be written to a single VCF for each variant class
    File? snv_partition_intervals
    File? indel_partition_intervals
    File? sv_partition_intervals

    String g2c_analysis_docker
  }

  Int min_large_indel_size = floor(min_sv_size / size_scalar)

  # Make indel/sv clustering intervals
  call MakeClusteringIntervals {
    input:
      genome_file = genome_file,
      nmask_bed = nmask_bed,
      min_nmask_size = ceil(2 * size_scalar * min_sv_size),
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Split SNVs, small indels, and large indels
  Array[Pair[File, File]] gatkhc_vcf_infos = zip(gatkhc_vcfs, gatkhc_vcf_idxs)
  scatter ( gatkhc_info in gatkhc_vcf_infos ) {
      call SplitGatkHcBySize {
        input:
          vcf = gatkhc_info.left,
          vcf_idx = gatkhc_info.right,
          large_indel_size = min_large_indel_size,
          g2c_analysis_docker = g2c_analysis_docker
      }
  }

  # Partition large indels across clustering intervals

  # Partition SVs across clustering intervals

  # Cluster large indels and SVs

  # Re-integrate small and large indels

  # Postprocess SNVs

  # Postprocess indels

  # Postprocess SVs
}


# Make intervals for indel/SV clustering based on N-masked reference regions
task MakeClusteringIntervals {
  input {
    File genome_file
    File nmask_bed
    Int min_nmask_size
    String g2c_analysis_docker
  }

  command <<<
    set -eu -o pipefail

    # Prep N-mask
    cut -f3 ~{nmask_bed} \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | bedtools merge -i - \
    | awk -v ms="$min_nmask_size" -v FS="\t" -v OFS="\t" \
      '{ if ($3-$2>=ms) print }' \
    | bgzip -c \
    > nmask.filtered.bed.gz

    # Subtract qualifying N-masked intervals from genome file
    awk -v FS="\t" -v OFS="\t" '{ print $1, 0, $2 }' \
    | bedtools subtract -a - -b nmask.filtered.bed.gz \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | awk -v OFS="\t" '{ print $0, "clustering_interval_"NR }' \
    | bgzip -c \
    > clustering.intervals.bed.gz
  >>>

  output {
    File intervals_bed = "clustering.intervals.bed.gz"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk 25 HDD"
    preemptible: 1
    maxRetries: 1
  }
}


# Split a single GATK-HC VCF by variant size
task SplitGatkHcBySize {
  input {
    File vcf
    File vcf_idx
    Int large_indel_size
    Float mem_gb = 1.75
    Int n_cpu = 1
    String g2c_analysis_docker
  }

  Int disk_gb = ceil(3 * size(vcf, "GB")) + 10
  String out_prefix = basename(vcf, ".vcf.gz")
  String snv_out_vcf = out_prefix + ".snv.vcf.gz"
  String si_out_vcf = out_prefix + ".small_indel.vcf.gz"
  String li_out_vcf = out_prefix + ".large_indel.vcf.gz"

  command <<<
    set -eu -o pipefail

    # Split input VCF
    /opt/pancan_germline_wgs/scripts/gatkhc_helpers/split_gatkhc_by_size.py \
      -i ~{vcf} \
      --large-indel-min-size ~{large_indel_size} \
      -o ~{out_prefix}

    # Tabix non-empty output VCFs and delete empty VCFs
    for vcf in "~{out_prefix}.*.vcf.gz"; do
      if [ $( bcftools index -n $vcf ) -gt 0 ]; then
        tabix -p vcf -f $vcf
      else
        rm $vcf
      fi
    done
  >>>

  output {
    File? snv_vcf = snv_out_vcf
    File? snv_vcf_idx = snv_out_vcf + ".tbi"

    File? small_indel_vcf = si_out_vcf
    File? small_indel_vcf_idx = si_out_vcf + ".tbi"

    File? large_indel_vcf = li_out_vcf
    File? large_indel_vcf_idx = li_out_vcf + ".tbi"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 1
    maxRetries: 1
  }
}

