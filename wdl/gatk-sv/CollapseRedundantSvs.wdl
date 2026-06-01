# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Collapse redundant SVs at the end of the GATK-SV pipeline
# Designed to be implemented after module 19 (FilterGenotypes)


version 1.0


workflow CollapseRedundantSvs {
  input {
    File vcf
    File vcf_idx

    Array[Float] recip_overlap = [0.975, 0.99, 1]
    Array[Int] bp_dist = [800, 200, 1]
    Array[Float] sample_overlap = [0.10, 0.05, 0]
    Array[Float] min_af = [0.01, 0.001, 0]
    Array[Int] min_ac = [100, 10, 1]
    Array[String] round_prefixes = ["loose", "strict", "identical"]

    String cluster_sv_types = "DEL,DUP,INS,INV,CPX"

    String g2c_analysis_docker
  }

  # ROUND 1
  call DefineClusters as DC1 {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      recip_overlap = recip_overlap[0],
      bp_dist = bp_dist[0],
      sample_overlap = sample_overlap[0],
      cluster_sv_types = cluster_sv_types,
      min_af = min_af[0],
      min_ac = min_ac[0],
      out_prefix = basename(vcf, "vcf.gz") + round_prefixes[0],
      g2c_analysis_docker = g2c_analysis_docker    
  }

  if ( DC1.n_clusters > 0 ) {
    call ResolveClusters as RC1 {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        clusters = DC1.clusters,
        out_prefix = basename(vcf, "vcf.gz") + round_prefixes[0],
        g2c_analysis_docker = g2c_analysis_docker
    }
  }
  File RC1_vcf = select_first([RC1.reclustered_vcf, vcf])
  File RC1_vcf_idx = select_first([RC1.reclustered_vcf_idx, vcf_idx])

  # ROUND 2
  call DefineClusters as DC2 {
    input:
      vcf = RC1_vcf,
      vcf_idx = RC1_vcf_idx,
      recip_overlap = recip_overlap[1],
      bp_dist = bp_dist[1],
      sample_overlap = sample_overlap[1],
      cluster_sv_types = cluster_sv_types,
      min_af = min_af[1],
      min_ac = min_ac[1],
      out_prefix = basename(RC1_vcf, "vcf.gz") + round_prefixes[1],
      g2c_analysis_docker = g2c_analysis_docker    
  }

  if ( DC2.n_clusters > 0 ) {
    call ResolveClusters as RC2 {
      input:
        vcf = RC1_vcf,
        vcf_idx = RC1_vcf_idx,
        clusters = DC2.clusters,
        out_prefix = basename(RC1_vcf, "vcf.gz") + round_prefixes[1],
        g2c_analysis_docker = g2c_analysis_docker
    }
  }
  File RC2_vcf = select_first([RC2.reclustered_vcf, RC1_vcf])
  File RC2_vcf_idx = select_first([RC2.reclustered_vcf_idx, RC1_vcf_idx])

  # ROUND 3
  call DefineClusters as DC3 {
    input:
      vcf = RC2_vcf,
      vcf_idx = RC2_vcf_idx,
      recip_overlap = recip_overlap[2],
      bp_dist = bp_dist[2],
      sample_overlap = sample_overlap[2],
      cluster_sv_types = cluster_sv_types,
      min_af = min_af[2],
      min_ac = min_ac[2],
      out_prefix = basename(RC2_vcf, "vcf.gz") + round_prefixes[2],
      g2c_analysis_docker = g2c_analysis_docker    
  }

  if ( DC3.n_clusters > 0 ) {
    call ResolveClusters as RC3 {
      input:
        vcf = RC2_vcf,
        vcf_idx = RC2_vcf_idx,
        clusters = DC3.clusters,
        out_prefix = basename(vcf, "vcf.gz") + round_prefixes[2],
        g2c_analysis_docker = g2c_analysis_docker
    }
  }

  output {
    File reclustered_vcf = select_first([RC3.reclustered_vcf, RC2_vcf])
    File reclustered_vcf_idx = select_first([RC3.reclustered_vcf_idx, RC2_vcf_idx])
    Array[File] cluster_maps = select_all([DC1.clusters, DC2.clusters, DC3.clusters])
  }
}


# Define variants to be clustered based on user-specified clustering parameters
task DefineClusters {
  input {
    File vcf
    File vcf_idx

    Float recip_overlap
    Int bp_dist
    Float sample_overlap
    String cluster_sv_types
    Float min_af
    Int min_ac

    String out_prefix
    
    Float mem_gb = 7.5
    Int n_cpu = 4

    String g2c_analysis_docker
  }

  Int disk_gb = ceil(2 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    # Assign variants to clusters using svtk
    echo "##source=gatksv" > add_header.vcf
    bcftools annotate -h add_header.vcf ~{vcf} \
    | bcftools +fill-tags -- -t AF,AC \
    | bcftools view --include "AF>=~{min_af} & AC>=~{min_ac}" -Oz -o input.vcf.gz
    tabix -f input.vcf.gz
    echo "input.vcf.gz" > input_vcf.list
    svtk vcfcluster \
      -d ~{bp_dist} \
      -f ~{recip_overlap} \
      -o ~{sample_overlap} \
      -p ~{out_prefix}.reclustered \
      -t "~{cluster_sv_types}" \
      --preserve-header \
      --skip-merge \
      input_vcf.list \
      - \
    | bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%INFO/CLUSTER\n' \
    > ~{out_prefix}.sv_cluster_assignments.bed

    # Identify clusters of two or more variants
    while read cidx; do
      awk -v cidx=$cidx -v OFS="\t" \
        '{ if ($5==cidx) print $1, $2, $3, $4 }' \
        ~{out_prefix}.sv_cluster_assignments.bed \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | bedtools merge -i - -d ~{bp_dist} -c 4 -o distinct
    done < <( cut -f5 ~{out_prefix}.sv_cluster_assignments.bed \
              | uniq -c | awk '{ if ($1>1) print $2 }' ) \
    | awk '{ if ($4~/,/) print }' \
    | sort -Vk1,1 -k2,2n -k3,3n \
    > ~{out_prefix}.sv_clusters.bed
  >>>

  output {
    File clusters = "~{out_prefix}.sv_clusters.bed"
    Int n_clusters = length(read_lines("~{out_prefix}.sv_clusters.bed"))
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "~{mem_gb} GB"
    cpu: n_cpu
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
    maxRetries: 1
  }
}


task ResolveClusters {
  input {
    File vcf
    File vcf_idx
    File clusters

    String out_prefix

    Float mem_gb = 7.5
    Int n_cpu = 4

    String g2c_analysis_docker
  }

  String vcf_out = "~{out_prefix}.reclustered.vcf.gz"

  Int disk_gb = ceil(10 * size(vcf, "GB")) + 10

  Float sort_mem = 2 * mem_gb / 3

  command <<<
    set -eu -o pipefail

    # Process each variant cluster and write to reclustered output file
    /opt/pancan_germline_wgs/scripts/gatksv_helpers/resolve_variant_clusters.py \
      --in-vcf ~{vcf} \
      --clusters ~{clusters} \
      --prefix "~{out_prefix}_reclustered" \
    | bcftools sort -m "~{sort_mem}G" \
    | bcftools +fill-tags -Oz -o ~{out_prefix}.resolved_clusters.vcf.gz -- -t AC,AN,AF
    tabix -p vcf -f ~{out_prefix}.resolved_clusters.vcf.gz

    # Remove cluster members from input VCF and merge remainder with new clustered records
    cut -f4 ~{clusters} | sed 's/,/\n/g' | sort | uniq > member.vids.list
    bcftools view --exclude 'ID=@member.vids.list' ~{vcf} -Oz -o ~{out_prefix}.remainder.vcf.gz
    tabix -p vcf -f ~{out_prefix}.remainder.vcf.gz
    bcftools concat \
      --allow-overlaps \
      ~{out_prefix}.remainder.vcf.gz \
      ~{out_prefix}.resolved_clusters.vcf.gz \
      -Oz -o ~{vcf_out} 
    tabix -p vcf -f ~{vcf_out}
  >>>

  output {
    File reclustered_vcf = vcf_out
    File reclustered_vcf_idx = "~{vcf_out}.tbi"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "~{mem_gb} GB"
    cpu: n_cpu
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
    maxRetries: 1
    bootDiskSizeGb: 30
  }
}
