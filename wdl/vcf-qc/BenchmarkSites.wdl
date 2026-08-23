# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Perform two-sided site-level benchmarking of a "source" callset versus a reference "target" callset

# Expects that both callsets have already been processed by CollectVcfQcMetrics.wdl


version 1.0


import "BenchmarkSitesSingle.wdl" as BenchSingle
import "QcTasks.wdl" as QcTasks


workflow BenchmarkSites {
  input {
    File? source_snv_bed
    File? source_indel_bed
    File? source_sv_bed
    String source_prefix

    File? target_snv_bed
    File? target_indel_bed
    File? target_sv_bed
    String target_prefix

    Array[File] eval_interval_beds
    Array[String]? eval_interval_bed_names
    File genome_file

    Int total_shards = 2500
    Float common_af_cutoff
    Int snv_error_run_min_variants = 10
    Int indel_error_run_min_variants = 10
    Int sv_error_run_min_variants = 10
    Int snv_error_run_min_bp = 20000
    Int indel_error_run_min_bp = 20000
    Int sv_error_run_min_bp = 50000

    Boolean make_full_id_maps = false

    Float snv_mem_scalar = 2.0
    Float indel_mem_scalar = 1.0
    Float shard_intervals_mem_gb = 3.75

    String bcftools_docker
    String g2c_analysis_docker
  }


  # Index SNV BEDs
  if ( defined(source_snv_bed) ){
    call QcTasks.MakeTabixIndex as IndexSourceSnvs {
      input:
        input_file = select_first([source_snv_bed]),
        file_type = "bed",
        docker = bcftools_docker
    }
  }
  if ( defined(target_snv_bed) ){
    call QcTasks.MakeTabixIndex as IndexTargetSnvs {
      input:
        input_file = select_first([target_snv_bed]),
        file_type = "bed",
        docker = bcftools_docker
    }
  }


  # Index indel BEDs
  if ( defined(source_indel_bed) ){
    call QcTasks.MakeTabixIndex as IndexSourceIndels {
      input:
        input_file = select_first([source_indel_bed]),
        file_type = "bed",
        docker = bcftools_docker
    }
  }
  if ( defined(target_indel_bed) ){
    call QcTasks.MakeTabixIndex as IndexTargetIndels {
      input:
        input_file = select_first([target_indel_bed]),
        file_type = "bed",
        docker = bcftools_docker
    }
  }


  # Index SV BEDs
  if ( defined(source_sv_bed) ){
    call QcTasks.MakeTabixIndex as IndexSourceSvs {
      input:
        input_file = select_first([source_sv_bed]),
        file_type = "bed",
        docker = bcftools_docker
    }
  }
  if ( defined(target_sv_bed) ){
    call QcTasks.MakeTabixIndex as IndexTargetSvs {
      input:
        input_file = select_first([target_sv_bed]),
        file_type = "bed",
        docker = bcftools_docker
    }
  }


  # Main benchmarking block
  Int n_eval_beds = length(eval_interval_beds)
  Int shards_per_eval_bed = floor(total_shards / n_eval_beds)

  scatter ( eval_idx in range(n_eval_beds) ){

    File eval_bed = eval_interval_beds[eval_idx]
    if ( defined(eval_interval_bed_names) ) {
      String eval_name_custom = select_first([eval_interval_bed_names])[eval_idx]
    }
    String eval_name = select_first([eval_name_custom, basename(eval_bed, ".bed.gz")])

    # Run evaluation for a single eval BED as a subworkflow
    call BenchSingle.BenchmarkSitesSingle as BenchmarkTask {
      input:
        eval_interval_bed = eval_bed,
        genome_file = genome_file,
        source_snv_bed = source_snv_bed,
        source_snv_bed_idx = IndexSourceSnvs.tbi,
        target_snv_bed = target_snv_bed,
        target_snv_bed_idx = IndexTargetSnvs.tbi,
        source_indel_bed = source_indel_bed,
        source_indel_bed_idx = IndexSourceIndels.tbi,
        target_indel_bed = target_indel_bed,
        target_indel_bed_idx = IndexTargetIndels.tbi,
        source_sv_bed = source_sv_bed,
        source_sv_bed_idx = IndexSourceSvs.tbi,
        target_sv_bed = target_sv_bed,
        target_sv_bed_idx = IndexTargetSvs.tbi,
        eval_prefix = eval_name,
        source_prefix = source_prefix,
        target_prefix = target_prefix,
        common_af_cutoff = common_af_cutoff,
        total_shards = shards_per_eval_bed,
        make_full_id_maps = make_full_id_maps,
        snv_mem_scalar = snv_mem_scalar,
        indel_mem_scalar = indel_mem_scalar,
        shard_intervals_mem_gb = shard_intervals_mem_gb,
        bcftools_docker = bcftools_docker,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }

  # Detect runs of false positive/negative SNVs
  Array[File] all_common_snv_ppv_beds = select_all(BenchmarkTask.common_snv_ppv_bed)
  if ( length(all_common_snv_ppv_beds) > 0 ) {
    call QcTasks.ConcatTextFiles as CollapseCommonSnvPpvBeds {
      input:
        shards = all_common_snv_ppv_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = "common_snv.ppv.bed.gz",
        docker = bcftools_docker
    }
    call DetectBadRuns as DetectSnvFPRuns {
      input:
        bed = CollapseCommonSnvPpvBeds.merged_file,
        output_prefix = "~{source_prefix}.~{target_prefix}.snv.false_positive_runs",
        min_variants = snv_error_run_min_variants,
        min_run_length = snv_error_run_min_bp,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }
  Array[File] all_common_snv_sens_beds = select_all(BenchmarkTask.common_snv_sens_bed)
  if ( length(all_common_snv_sens_beds) > 0 ) {
    call QcTasks.ConcatTextFiles as CollapseCommonSnvSensBeds {
      input:
        shards = all_common_snv_sens_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = "common_snv.sens.bed.gz",
        docker = bcftools_docker
    }
    call DetectBadRuns as DetectSnvFNRuns {
      input:
        bed = CollapseCommonSnvSensBeds.merged_file,
        output_prefix = "~{source_prefix}.~{target_prefix}.snv.false_negative_runs",
        min_variants = snv_error_run_min_variants,
        min_run_length = snv_error_run_min_bp,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }


  # Detect runs of false positive/negative indels
  Array[File] all_common_indel_ppv_beds = select_all(BenchmarkTask.common_indel_ppv_bed)
  if ( length(all_common_indel_ppv_beds) > 0 ) {
    call QcTasks.ConcatTextFiles as CollapseCommonIndelPpvBeds {
      input:
        shards = all_common_indel_ppv_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = "common_indel.ppv.bed.gz",
        docker = bcftools_docker
    }
    call DetectBadRuns as DetectIndelFPRuns {
      input:
        bed = CollapseCommonIndelPpvBeds.merged_file,
        output_prefix = "~{source_prefix}.~{target_prefix}.indel.false_positive_runs",
        min_variants = indel_error_run_min_variants,
        min_run_length = indel_error_run_min_bp,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }
  Array[File] all_common_indel_sens_beds = select_all(BenchmarkTask.common_indel_sens_bed)
  if ( length(all_common_indel_sens_beds) > 0 ) {
    call QcTasks.ConcatTextFiles as CollapseCommonIndelSensBeds {
      input:
        shards = all_common_indel_sens_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = "common_indel.sens.bed.gz",
        docker = bcftools_docker
    }
    call DetectBadRuns as DetectIndelFNRuns {
      input:
        bed = CollapseCommonIndelSensBeds.merged_file,
        output_prefix = "~{source_prefix}.~{target_prefix}.indel.false_negative_runs",
        min_variants = indel_error_run_min_variants,
        min_run_length = indel_error_run_min_bp,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }


  # Detect runs of false positive/negative SVs
  Array[File] all_common_sv_ppv_beds = select_all(BenchmarkTask.common_sv_ppv_bed)
  if ( length(all_common_sv_ppv_beds) > 0 ) {
    call QcTasks.ConcatTextFiles as CollapseCommonSvPpvBeds {
      input:
        shards = all_common_sv_ppv_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = "common_sv.ppv.bed.gz",
        docker = bcftools_docker
    }
    call DetectBadRuns as DetectSvFPRuns {
      input:
        bed = CollapseCommonSvPpvBeds.merged_file,
        output_prefix = "~{source_prefix}.~{target_prefix}.sv.false_positive_runs",
        min_variants = sv_error_run_min_variants,
        min_run_length = sv_error_run_min_bp,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }
  Array[File] all_common_sv_sens_beds = select_all(BenchmarkTask.common_sv_sens_bed)
  if ( length(all_common_sv_sens_beds) > 0 ) {
    call QcTasks.ConcatTextFiles as CollapseCommonSvSensBeds {
      input:
        shards = all_common_sv_sens_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = "common_sv.sens.bed.gz",
        docker = bcftools_docker
    }
    call DetectBadRuns as DetectSvFNRuns {
      input:
        bed = CollapseCommonSvSensBeds.merged_file,
        output_prefix = "~{source_prefix}.~{target_prefix}.sv.false_negative_runs",
        min_variants = sv_error_run_min_variants,
        min_run_length = sv_error_run_min_bp,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }


  output {
    Array[File] common_snv_ppv_beds = select_all(BenchmarkTask.common_snv_ppv_bed)
    Array[File] common_snv_sens_beds = select_all(BenchmarkTask.common_indel_ppv_bed)
    Array[File] common_indel_ppv_beds = select_all(BenchmarkTask.common_indel_ppv_bed)
    Array[File] common_indel_sens_beds = select_all(BenchmarkTask.common_indel_sens_bed)
    Array[File] common_sv_ppv_beds = select_all(BenchmarkTask.common_sv_ppv_bed)
    Array[File] common_sv_sens_beds = select_all(BenchmarkTask.common_sv_sens_bed)
    Array[File] ppv_by_sizes = BenchmarkTask.ppv_by_size
    Array[File] sensitivity_by_sizes = BenchmarkTask.sensitivity_by_size
    Array[File] ppv_by_freqs = BenchmarkTask.ppv_by_freq
    Array[File] sensitivity_by_freqs = BenchmarkTask.sensitivity_by_freq
    Array[File?] ppv_variant_id_maps = BenchmarkTask.ppv_variant_id_map
    Array[File?] sensitivity_variant_id_maps = BenchmarkTask.sensitivity_variant_id_map
    Array[File] fp_runs = select_all([DetectSnvFPRuns.runs_bed, DetectIndelFPRuns.runs_bed, DetectSvFPRuns.runs_bed])
    Array[File] fn_runs = select_all([DetectSnvFNRuns.runs_bed, DetectIndelFNRuns.runs_bed, DetectSvFNRuns.runs_bed])
  }
}


task DetectBadRuns {
  input {
    File bed
    String output_prefix

    Int min_variants = 10
    Int min_run_length = 50000

    String g2c_analysis_docker
    String gcp_machine_type = "n2d-standard-2"
  }

  Int disk_gb = ceil(2 * size(bed, "GB")) + 10
  String out_fname = "~{output_prefix}.bed.gz"

  command <<<
    set -eu -o pipefail

    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/find_problematic_intervals.py \
      --min-vars ~{min_variants} \
      --min-bp ~{min_run_length} \
      ~{bed} \
      ~{out_fname}
  >>>

  output {
    File runs_bed = out_fname
  }

  runtime {
    docker: g2c_analysis_docker
    predefinedMachineType: gcp_machine_type
    memory: "3.7 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 2
    max_retries: 1
  }
}

