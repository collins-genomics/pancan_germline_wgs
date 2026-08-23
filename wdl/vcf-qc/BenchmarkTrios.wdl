# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Perform Mendelian violation analysis of parent-child trios


version 1.0


import "QcTasks.wdl" as QcTasks


workflow BenchmarkTrios {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    File sites_bed
    Boolean concat_vcfs = false

    File trios_fam
    File trios_samples_list

    Array[File] eval_interval_beds
    Array[String]? eval_interval_bed_names
    File genome_file

    String output_prefix
    Float common_af_cutoff = 0.01

    Float benchmarking_mem_gb = 3.75
    Int benchmarking_n_cpu = 2

    String bcftools_docker
    String g2c_analysis_docker
  }

  # Index sites BED
  call QcTasks.MakeTabixIndex as IndexSites {
    input:
      input_file = sites_bed,
      file_type = "bed",
      docker = bcftools_docker
  }

  # Concatenate VCF shards, if optioned
  if ( concat_vcfs ) {
    call QcTasks.ConcatVcfs as ConcatShards {
      input:
        vcfs = vcfs,
        vcf_idxs = vcf_idxs,
        out_prefix = output_prefix + ".merged_shards_for_trio_bench",
        bcftools_concat_options = "--allow-overlaps",
        bcftools_docker = bcftools_docker
    }
  }
  Array[File] vcfs_use = if concat_vcfs then select_all([ConcatShards.merged_vcf]) else vcfs
  Array[File] vcf_idxs_use = if concat_vcfs then select_all([ConcatShards.merged_vcf_idx]) else vcf_idxs

  # Keep benchmarking results separate for each evaluation interval set
  Int n_eval_intervals = length(eval_interval_beds)
  scatter ( i in range(n_eval_intervals) ) {

    String ei_prefix = if defined(eval_interval_bed_names)
                       then flatten(select_all([eval_interval_bed_names]))[i]
                       else "eval_interval_~{i}"

    # Subset sites to those within eval intervals
    call QcTasks.PrepSites as SliceSites {
      input:
        beds = [sites_bed],
        bed_idxs = [IndexSites.tbi],
        eval_interval_bed = eval_interval_beds[i],
        prefix = ei_out_prefix,
        g2c_analysis_docker = g2c_analysis_docker
    }

    # Benchmark trios in each VCF
    Array[Pair[File, File]] vcf_infos = zip(vcfs_use, vcf_idxs_use)
    scatter ( vcf_info in vcf_infos ) {
      call BenchmarkTrios {
        input:
          vcf = vcf_info.left,
          vcf_idx = vcf_info.right,
          eligible_sites_bed = SliceSites.query_bed,
          trios_fam = trios_fam,
          trios_samples_list = trios_samples_list,
          common_af_cutoff = common_af_cutoff,
          mem_gb = benchmarking_mem_gb,
          n_cpu = benchmarking_n_cpu,
          output_prefix = ei_prefix + "." + basename(vcf_info.left, ".vcf.gz"),
          g2c_analysis_docker = g2c_analysis_docker
      }
    }

    # Collapse results across all VCFs within each evaluation interval set
    String ei_out_prefix = "~{ei_prefix}.~{output_prefix}"
    call QcTasks.SumCompressedDistribs {
      input:
        distrib_tsvs = BenchmarkTrios.benchmarking_distrib,
        n_key_columns = 4,
        out_prefix = ei_out_prefix + ".mendelian_distribution",
        g2c_analysis_docker = g2c_analysis_docker
    }
  }

  output {
    Array[File] trio_benchmarking_distribs = SumCompressedDistribs.merged_distrib
  }
}


task BenchmarkTrios {
  input {
    File vcf
    File vcf_idx
    File eligible_sites_bed
    
    File trios_fam
    File trios_samples_list

    Float common_af_cutoff = 0.01
    String output_prefix
    
    String gcp_machine_type = "n2d-standard-2"
    Float mem_gb
    Int n_cpu
    String g2c_analysis_docker
  }

  String out_fname = "~{output_prefix}.mendelian_violations.distribs.tsv.gz"

  String machine_type = if defined(gcp_machine_type) then gcp_machine_type else ""
  Int disk_gb = ceil(3 * size([vcf, eligible_sites_bed], "GB")) + 20

  command <<<
    set -eu -o pipefail

    # Gather list of eligible variant IDs
    zcat ~{eligible_sites_bed} | fgrep -v "#" \
    | cut -f4 | sort -V | uniq \
    > elig_vids.list

    # Subset VCF to samples from trios and benchmark Mendelian violations
    # Don't do this in one single pipe to avoid VM timeouts for very large VCFs
    bcftools view --no-update --force-samples \
      --samples-file ~{trios_samples_list} \
      --include 'ID==@elig_vids.list' \
      ~{vcf} \
    | bcftools view --no-update \
      --include 'GT="alt" | FILTER="MULTIALLELIC"' \
      -Oz -o input.vcf.gz
    tabix -p vcf -f input.vcf.gz
    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/gather_mendelian_violations.py \
      --trios-fam ~{trios_fam} \
      --eligible-vids elig_vids.list \
      --common-af ~{common_af_cutoff} \
      --summary-out stdout \
      input.vcf.gz \
    | gzip -c > ~{out_fname}
  >>>

  output {
    File benchmarking_distrib = out_fname
  }

  runtime {
    docker: g2c_analysis_docker
    predefinedMachineType: machine_type
    memory: "~{mem_gb} GB"
    cpu: n_cpu
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}

