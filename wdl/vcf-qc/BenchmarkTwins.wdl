# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Perform sample-level benchmarking of technical replicates/twins


version 1.0


import "QcTasks.wdl" as QcTasks


workflow BenchmarkTwins {
  input {
    File all_sites_bed
    File gt_tarball
    File twins_tsv
    String output_prefix

    Array[File] eval_interval_beds
    Array[String]? eval_interval_bed_names
    File genome_file

    Float common_af_cutoff = 0.01

    String bcftools_docker
    String g2c_analysis_docker
  }

  # Index sites BED
  call QcTasks.MakeTabixIndex as IndexSites {
    input:
      input_file = all_sites_bed,
      file_type = "bed",
      docker = bcftools_docker
  }

  # Keep benchmarking results separate for each evaluation interval set
  Int n_eval_intervals = length(eval_interval_beds)
  scatter ( i in range(n_eval_intervals) ) {

    String ei_prefix = if defined(eval_interval_bed_names)
                       then flatten(select_all([eval_interval_bed_names]))[i]
                       else "eval_interval_~{i}"
    String ei_out_prefix = "~{ei_prefix}.~{output_prefix}"

    # Subset sites to those within eval intervals
    call QcTasks.PrepSites as SliceSites {
      input:
        beds = [all_sites_bed],
        bed_idxs = [IndexSites.tbi],
        eval_interval_bed = eval_interval_beds[i],
        prefix = ei_out_prefix,
        g2c_analysis_docker = g2c_analysis_docker
    }

    # Make a fake variant ID map for all variants in these eval intervals
    call MakeVidMap {
      input:
        sites_bed = SliceSites.query_bed,
        prefix = output_prefix
    }

    # Benchmark sample 1 vs. sample 2
    call QcTasks.BenchmarkGenotypes as BenchGt1v2 {
      input:
        variant_id_map = MakeVidMap.vid_map_tsv,
        sample_id_map = twins_tsv,
        invert_sample_map = false,
        output_prefix = "~{ei_out_prefix}.1v2",
        source_gt_tarball = gt_tarball,
        target_gt_tarball = gt_tarball,
        source_site_metrics = SliceSites.query_bed,
        report_by_gt = false,
        common_af_cutoff = common_af_cutoff,
        g2c_analysis_docker = g2c_analysis_docker
    }

    # Benchmark sample 2 vs. sample 1
    call QcTasks.BenchmarkGenotypes as BenchGt2v1 {
      input:
        variant_id_map = MakeVidMap.vid_map_tsv,
        sample_id_map = twins_tsv,
        invert_sample_map = true,
        output_prefix = "~{ei_out_prefix}.2v1",
        source_gt_tarball = gt_tarball,
        target_gt_tarball = gt_tarball,
        source_site_metrics = SliceSites.query_bed,
        report_by_gt = false,
        common_af_cutoff = common_af_cutoff,
        g2c_analysis_docker = g2c_analysis_docker
    }

    # Collapse two benchmarking directions
    call QcTasks.SumCompressedDistribs {
      input:
        distrib_tsvs = [BenchGt1v2.gt_bench_distrib, BenchGt2v1.gt_bench_distrib],
        out_prefix = ei_out_prefix + ".gt_comparison.distrib",
        n_key_columns = 5,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }

  output {
    Array[File] gt_benchmarking_distribs = SumCompressedDistribs.merged_distrib
  }
}


task MakeVidMap {
  input {
    File sites_bed
    String prefix
    String linux_docker = "ubuntu:plucky-20251001"
  }

  String vid_map_fname = "~{prefix}.vid_map.tsv.gz"
  Int disk_gb = ceil(3 * size(sites_bed, "GB")) + 10

  command <<<
    set -eu -o -pipefail

    zcat ~{sites_bed} \
    | fgrep -v "#" \
    | awk -v OFS="\t" '{ print $4, $4 }' \
    | gzip -c \
    > ~{vid_map_fname}
  >>>

  output {
    File vid_map_tsv = "~{vid_map_fname}"
  }

  runtime {
    docker: linux_docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }  
}
