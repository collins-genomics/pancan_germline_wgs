# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Perform two-sided sample-level benchmarking of a "source" callset versus a reference "target" callset

# Expects that both callsets comply to formatting expected from the G2C QC workflow


version 1.0


import "BenchmarkSamplesSingle.wdl" as BenchSingle
import "QcTasks.wdl" as QcTasks


workflow BenchmarkSamples {
  input {
    File? source_all_sites_bed
    Array[File?] source_snv_beds
    Array[File?] source_indel_beds
    Array[File?] source_sv_beds
    File source_gt_tarball
    String source_prefix

    Array[File] target_vcfs
    Array[File] target_vcf_idxs
    String target_prefix

    Array[File] id_map_tsvs

    Array[File] eval_interval_beds
    Array[String]? eval_interval_bed_names
    File genome_file

    Int total_shards = 2500
    Int min_samples_per_shard = 10
    Float common_af_cutoff = 0.01

    Float snv_mem_scalar = 2.0
    Float indel_mem_scalar = 1.0

    String bcftools_docker
    String g2c_analysis_docker
  }

  Int n_targets = length(target_vcfs)
  Int n_shards_per_vcf = ceil(total_shards / n_targets)
  Int n_eval_intervals = length(eval_interval_beds)

  Array[File] source_snv_beds_nonull = select_all(source_snv_beds)
  Array[File] source_indel_beds_nonull = select_all(source_indel_beds)
  Array[File] source_sv_beds_nonull = select_all(source_sv_beds)

  call QcTasks.MakeEmptyBenchBed {
    input:
      docker = bcftools_docker
  }

  # Collapse source SNV sites
  if ( length(source_snv_beds_nonull) > 0 ) {
    call QcTasks.ConcatTextFiles as CollapseSourceSnvs {
      input:
        shards = source_snv_beds_nonull,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = source_prefix + ".all_snvs.bed.gz",
        docker = bcftools_docker
    }
  }
  File source_snv_bed = select_first([CollapseSourceSnvs.merged_file, MakeEmptyBenchBed.empty_bed])

  # Collapse source indel sites
  if ( length(source_indel_beds_nonull) > 0 ) {
    call QcTasks.ConcatTextFiles as CollapseSourceIndels {
      input:
        shards = source_indel_beds_nonull,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = source_prefix + ".all_indels.bed.gz",
        docker = bcftools_docker
    }
  }
  File source_indel_bed = select_first([CollapseSourceIndels.merged_file, MakeEmptyBenchBed.empty_bed])

  # Collapse source SV sites
  if ( length(source_sv_beds_nonull) > 0 ) {
    call QcTasks.ConcatTextFiles as CollapseSourceSvs {
      input:
        shards = source_sv_beds_nonull,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = source_prefix + ".all_svs.bed.gz",
        docker = bcftools_docker
    }
  }
  File source_sv_bed = select_first([CollapseSourceSvs.merged_file, MakeEmptyBenchBed.empty_bed])

  # Collapse all source site metrics unless already provided
  if ( !defined(source_all_sites_bed) ) {
    call QcTasks.ConcatTextFiles as CollapseAllSourceVars {
      input:
        shards = select_all([source_snv_bed, source_indel_bed, source_sv_bed]),
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = source_prefix + ".all_sites.bed.gz",
        docker = bcftools_docker
    }
  }
  File source_all_sites = select_first([source_all_sites_bed, CollapseAllSourceVars.merged_file])
  

  # Determine sample IDs present in GT tarball
  call GetSampleIdsFromGtTarball {
    input:
      gt_tarball = source_gt_tarball,
      outfile_name = "source_sample_ids.list"
  }

  # Benchmark each target VCF in parallel
  scatter ( i in range(n_targets) ) {
    call BenchSingle.BenchmarkSamplesSingle as BenchmarkTask {
      input:
        source_all_sites_bed = source_all_sites,
        source_snv_bed = source_snv_bed,
        source_indel_bed = source_indel_bed,
        source_sv_bed = source_sv_bed,
        source_gt_tarball = source_gt_tarball,
        source_sample_id_list = GetSampleIdsFromGtTarball.sample_id_list,
        source_prefix = source_prefix,
        target_vcf = target_vcfs[i],
        target_vcf_idx = target_vcf_idxs[i],
        target_prefix = target_prefix,
        id_map_tsv = id_map_tsvs[i],
        eval_interval_beds = eval_interval_beds,
        eval_interval_bed_names = eval_interval_bed_names,
        genome_file = genome_file,
        total_shards = n_shards_per_vcf,
        min_samples_per_shard = min_samples_per_shard,
        common_af_cutoff = common_af_cutoff,
        snv_mem_scalar = snv_mem_scalar,
        indel_mem_scalar = indel_mem_scalar,
        bcftools_docker = bcftools_docker,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }

  # BenchmarkTask is scattered per input VCF, and each BenchmarkTask call
  # outputs an array with one element per eval interval set. Since we want 
  # to collapse across input VCFs to produce one composite file per eval 
  # interval, we need to transpose the results from BenchmarkTask
  Array[Array[File]] ppv_distribs_t = transpose(BenchmarkTask.ppv_distribs)
  Array[Array[File]] sens_distribs_t = transpose(BenchmarkTask.sensitivity_distribs)

  # Pool benchmarking results across VCFs per evalulation interval set
  scatter ( k in range(n_eval_intervals) ) {

    String ei_prefix = if defined(eval_interval_bed_names)
                       then flatten(select_all([eval_interval_bed_names]))[k]
                       else "eval_interval_~{k}"
    String ppv_prefix = "~{ei_prefix}.~{source_prefix}.~{target_prefix}"
    String sens_prefix = "~{ei_prefix}.~{target_prefix}.~{source_prefix}"

    # Collapse PPV benchmarking across all target VCFs
    call QcTasks.SumCompressedDistribs as SumPpvDistrib {
      input:
        distrib_tsvs = ppv_distribs_t[k],
        out_prefix = ppv_prefix + ".gt_comparison.distrib",
        n_key_columns = 5,
        g2c_analysis_docker = g2c_analysis_docker
    }

    # Collapse sensitivity benchmarking across all target VCFs
    call QcTasks.SumCompressedDistribs as SumSensDistrib {
      input:
        distrib_tsvs = sens_distribs_t[k],
        out_prefix = sens_prefix + ".gt_comparison.distrib",
        n_key_columns = 5,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }

  output {
    # One file per evaluation interval set in each output array
    Array[File] ppv_distribs = SumPpvDistrib.merged_distrib
    Array[File] sensitivity_distribs = SumSensDistrib.merged_distrib
  }
}


# Get the list of sample IDs present in a GT tarball generated by QcTasks.ConcatGenotypeTsvs
task GetSampleIdsFromGtTarball {
  input {
    File gt_tarball
    String outfile_name = "ids_in_tarball.list"
    String linux_docker = "ubuntu:plucky-20251001"
  }

  Int disk_gb = ceil(4 * size(gt_tarball, "GB")) + 10

  command <<<
    set -eu -o pipefail

    mkdir gt_out

    tar -xzvf ~{gt_tarball} -C gt_out/

    find gt_out/ -name "*.gt.tsv.gz" \
    | xargs -I {} basename {} \
    | sed 's/\.gt\.tsv\.gz//g' \
    | sort \
    > ~{outfile_name} || true
  >>>

  output {
    File sample_id_list = outfile_name
  }

  runtime {
    docker: linux_docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}

