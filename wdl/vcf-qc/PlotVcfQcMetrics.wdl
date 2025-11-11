# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Plot quality metrics for assessing a GATK joint-genotyped VCF (SNVs/indels or SVs)
# See CollectVcfQcMetrics.wdl for the generation of inputs for this workflow


version 1.0


import "PrepSiteBenchDataToPlot.wdl" as PSB
import "QcTasks.wdl" as QcTasks


workflow PlotVcfQcMetrics {
  input {
    Array[File] size_distribution_tsvs
    Array[File] af_distribution_tsvs
    Array[File] size_vs_af_distribution_tsvs
    Array[File] sample_genotype_distribution_tsvs
    Array[File]? peak_ld_stat_tsvs

    Array[File]? ref_size_distribution_tsvs
    Array[File]? ref_af_distribution_tsvs
    String ref_cohort_prefix = "ref_dataset"
    String ref_cohort_plot_title = "Ref. data"

    Array[File]? all_sv_beds
    Array[File]? common_snv_beds
    Array[File]? common_indel_beds
    Array[File]? common_sv_beds

    File? sample_ancestry_labels
    File? sample_phenotype_labels
    File? sample_subset_list
    Float common_af_cutoff = 0.01
    Array[String?] benchmark_interval_names = []

    # Expected organization of site benchmarking inputs:
    # Outer array: one entry per benchmarking dataset
    # Middle arrays: one entry per evaluation interval set
    # Inner arrays: one or more files, such as one file per chromosome 
    # (or just one file if already collapsed across all chromosomes)
    Array[Array[Array[File?]]] site_benchmark_common_snv_ppv_beds = [[[]]]
    Array[Array[Array[File?]]] site_benchmark_common_snv_sens_beds = [[[]]]
    Array[Array[Array[File?]]] site_benchmark_common_indel_ppv_beds = [[[]]]
    Array[Array[Array[File?]]] site_benchmark_common_indel_sens_beds = [[[]]]
    Array[Array[Array[File?]]] site_benchmark_common_sv_ppv_beds = [[[]]]
    Array[Array[Array[File?]]] site_benchmark_common_sv_sens_beds = [[[]]]
    Array[Array[Array[File?]]] site_benchmark_ppv_by_freqs = [[[]]]
    Array[Array[Array[File?]]] site_benchmark_sensitivity_by_freqs = [[[]]]
    Array[String?] site_benchmark_dataset_prefixes = []
    Array[String?] site_benchmark_dataset_titles = []
    # For platforms like Terra, it is often more practical to provide inputs 
    # with middle & inner arrays flipped (i.e., middle array = per chromosome,
    # inner array = per eval set). If optioned, the below argument will 
    # automatically transpose the inner & middle arrays of all of the
    # site_benchmark_* inputs above to revert to the expected order
    Boolean transpose_site_benchmarking_nested_arrays = false

    # Expected organization of external sample-level genotype benchmarking inputs:
    # Outer array: one entry per benchmarking dataset
    # Middle arrays: one entry per evaluation interval set
    # Inner arrays: one or more files, such as one file per chromosome 
    # (or just one file if already collapsed across all chromosomes)
    Array[Array[Array[File?]]] sample_benchmark_ppv_distribs = [[[]]]
    Array[Array[Array[File?]]] sample_benchmark_sensitivity_distribs = [[[]]]
    Array[String?] sample_benchmark_dataset_prefixes = []
    Array[String?] sample_benchmark_dataset_titles = []
    # We also provide an option to invert the order of middle/inner arrays here
    # See commentary above for site_benchmark_* inputs
    Boolean transpose_sample_benchmarking_nested_arrays = false

    # Other benchmarking data; one outer array per evaluation interval
    # and each inner array can be one or more files (e.g., one per chromosome)
    Array[Array[File?]] twin_genotype_benchmark_distribs = [[]]
    Array[Array[File?]] trio_mendelian_violation_distribs = [[]]

    # Optionally, you can provide the all_stats_tsv output from a prior run for comparison
    File? previous_stats
    File? custom_qc_target_metrics
    File? custom_plotting_constants

    String output_prefix

    String bcftools_docker
    String g2c_analysis_docker
  }

  Int n_bench_intervals = length(select_first([benchmark_interval_names]))

  # Postprocess site benchmarking inputs
  Int n_sb_snv_beds = length(flatten(flatten(site_benchmark_common_snv_ppv_beds)))
  Int n_sb_indel_beds = length(flatten(flatten(site_benchmark_common_indel_ppv_beds)))
  Int n_sb_sv_beds = length(flatten(flatten(site_benchmark_common_sv_ppv_beds)))
  Int n_sb_ppv_by_afs = length(flatten(flatten(site_benchmark_ppv_by_freqs)))
  Int n_sb_sens_by_afs = length(flatten(flatten(site_benchmark_sensitivity_by_freqs)))
  Int n_sb_datasets = length(select_first([site_benchmark_dataset_prefixes]))
  Boolean has_site_benchmarking = ( (n_sb_snv_beds > 0 
                                     || n_sb_indel_beds > 0
                                     || n_sb_sv_beds > 0
                                     || n_sb_ppv_by_afs > 0
                                     || n_sb_sens_by_afs > 0) 
                                    && n_sb_datasets > 0 
                                    && n_bench_intervals > 0)

  # Postprocess external sample-level benchmarking inputs
  Int n_gb_ppv_tsvs = length(flatten(flatten(sample_benchmark_ppv_distribs)))
  Int n_gb_sens_tsvs = length(flatten(flatten(sample_benchmark_sensitivity_distribs)))
  Int n_gb_datasets = length(select_first([sample_benchmark_dataset_prefixes]))
  Boolean has_sample_benchmarking = ( (n_gb_ppv_tsvs > 0 || n_gb_sens_tsvs > 0) 
                                      && n_gb_datasets > 0 
                                      && n_bench_intervals > 0)

  # Postprocess other sample-level benchmarking (twins, trios)
  Int n_twin_bench_distribs = length(flatten(twin_genotype_benchmark_distribs))
  Boolean has_twin_benchmarking = (n_twin_bench_distribs > 0 && n_bench_intervals > 0)
  Int n_trio_bench_distribs = length(flatten(trio_mendelian_violation_distribs))
  Boolean has_trio_benchmarking = (n_trio_bench_distribs > 0 && n_bench_intervals > 0)


  ######################
  ### DATA PREPROCESSING
  ######################

  # If necessary, collapse size distribution into a single file
  if ( length(size_distribution_tsvs) > 1 ) {
      call QcTasks.SumCompressedDistribs as SumSizeDistribs {
        input:
          distrib_tsvs = size_distribution_tsvs,
          out_prefix = output_prefix + ".size_distribution",
          g2c_analysis_docker = g2c_analysis_docker
    }
  }
  File size_distrib = select_first([SumSizeDistribs.merged_distrib, size_distribution_tsvs[0]])

  # If necessary, collapse AF distribution into a single file
  if ( length(af_distribution_tsvs) > 1 ) {
      call QcTasks.SumCompressedDistribs as SumAfDistribs {
        input:
          distrib_tsvs = af_distribution_tsvs,
          out_prefix = output_prefix + ".af_distribution",
          g2c_analysis_docker = g2c_analysis_docker
    }
  }
  File af_distrib = select_first([SumAfDistribs.merged_distrib, af_distribution_tsvs[0]])

  # If necessary, collapse size-vs-AF distribution into a single file
  if ( length(size_vs_af_distribution_tsvs) > 1 ) {
      call QcTasks.SumCompressedDistribs as SumJointDistribs {
        input:
          distrib_tsvs = size_vs_af_distribution_tsvs,
          n_key_columns = 3,
          out_prefix = output_prefix + ".size_vs_af_distribution",
          g2c_analysis_docker = g2c_analysis_docker
    }
  }
  File joint_distrib = select_first([SumJointDistribs.merged_distrib, size_vs_af_distribution_tsvs[0]])

  # If necessary, collapse reference size distribution into a single file
  if ( defined(ref_size_distribution_tsvs) ) {
    Array[File] rsdt_array = select_all(select_first([ref_size_distribution_tsvs]))
    if ( length(rsdt_array) > 1 ) {
      call QcTasks.SumCompressedDistribs as SumRefSizeDistribs {
        input:
          distrib_tsvs = rsdt_array,
          out_prefix = ref_cohort_prefix + ".size_distribution",
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
    File ref_size_distrib = select_first([SumRefSizeDistribs.merged_distrib, rsdt_array[0]])
  }

  # If necessary, collapse reference AF distribution into a single file
  if ( defined(ref_af_distribution_tsvs) ) {
    Array[File] radt_array = select_all(select_first([ref_af_distribution_tsvs]))
    if ( length(radt_array) > 1 ) {
      call QcTasks.SumCompressedDistribs as SumRefAfDistribs {
        input:
          distrib_tsvs = radt_array,
          out_prefix = ref_cohort_prefix + ".af_distribution",
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
    File ref_af_distrib = select_first([SumRefAfDistribs.merged_distrib, radt_array[0]])
  }

  # If necessary, collapse all SV BEDs
  if ( defined(all_sv_beds) ) {
    if ( length(select_first(select_all([all_sv_beds]))) > 1 ) {
      call QcTasks.ConcatTextFiles as CollapseAllSvs {
        input:
          shards = select_first([all_sv_beds]),
          concat_command = "zcat",
          sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
          compression_command = "bgzip -c",
          input_has_header = true,
          output_filename = output_prefix + ".all_svs.bed.gz",
          docker = bcftools_docker
      }
    }
    File all_svs_bed = select_first(select_all([CollapseAllSvs.merged_file, 
                                                select_first([all_sv_beds])[0]]))
  }

  # If necessary, collapse common SNV BEDs
  if ( defined(common_snv_beds) ) {
    if ( length(select_first(select_all([common_snv_beds]))) > 1 ) {
      call QcTasks.ConcatTextFiles as CollapseCommonSnvs {
        input:
          shards = select_first([common_snv_beds]),
          concat_command = "zcat",
          sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
          compression_command = "bgzip -c",
          input_has_header = true,
          output_filename = output_prefix + ".common_snvs.bed.gz",
          docker = bcftools_docker
        }
    }
    File common_snvs_bed = select_first(select_all([CollapseCommonSnvs.merged_file, 
                                                    select_first([common_snv_beds])[0]]))
  }

  # If necessary, collapse common indel BEDs
  if ( defined(common_indel_beds) ) {
    if ( length(select_first(select_all([common_indel_beds]))) > 1 ) {
      call QcTasks.ConcatTextFiles as CollapseCommonIndels {
        input:
          shards = select_first([common_indel_beds]),
          concat_command = "zcat",
          sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
          compression_command = "bgzip -c",
          input_has_header = true,
          output_filename = output_prefix + ".common_indels.bed.gz",
          docker = bcftools_docker
        }
    }
    File common_indels_bed = select_first(select_all([CollapseCommonIndels.merged_file, 
                                                      select_first([common_indel_beds])[0]]))
  }

  # If necessary, collapse common SV BEDs
  if ( defined(common_sv_beds) ) {
    if ( length(select_first(select_all([common_sv_beds]))) > 1 ) {
      call QcTasks.ConcatTextFiles as CollapseCommonSvs {
        input:
          shards = select_first([common_sv_beds]),
          concat_command = "zcat",
          sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
          compression_command = "bgzip -c",
          input_has_header = true,
          output_filename = output_prefix + ".common_svs.bed.gz",
          docker = bcftools_docker
        }
    }
    File common_svs_bed = select_first(select_all([CollapseCommonSvs.merged_file, 
                                                   select_first([common_sv_beds])[0]]))
  }

  # If necessary, collapse sample genotype distributions
  if ( length(sample_genotype_distribution_tsvs) > 1 ) {
      call QcTasks.SumCompressedDistribs as SumGtDistribs {
        input:
          distrib_tsvs = sample_genotype_distribution_tsvs,
          out_prefix = output_prefix + ".genotype_distribution",
          n_key_columns = 4,
          g2c_analysis_docker = g2c_analysis_docker
    }
  }
  File gt_distrib = select_first([SumGtDistribs.merged_distrib, sample_genotype_distribution_tsvs[0]])

  # If necessary, collapse peak LD stats
  if ( defined(peak_ld_stat_tsvs) ) {
    Array[File] peak_ld_stat_tsv_use = select_first([peak_ld_stat_tsvs])
    if ( length(peak_ld_stat_tsv_use) > 1 ) {
      call QcTasks.ConcatTextFiles as CollapseLdStats {
        input:
          shards = peak_ld_stat_tsv_use,
          concat_command = "zcat",
          sort_command = "sort -Vk1,1 -k2,2V",
          compression_command = "gzip -c",
          input_has_header = true,
          output_filename = output_prefix + ".peak_ld_stats.tsv.gz",
          docker = bcftools_docker
      }
    }
    File ld_stats_tsv = select_first([CollapseLdStats.merged_file, peak_ld_stat_tsv_use[0]])
  }

  # Preprocess site benchmarking, if provided
  if (has_site_benchmarking) {
    scatter ( site_bench_di in range(n_sb_datasets) ) {

      String bd_name = select_first([site_benchmark_dataset_prefixes[site_bench_di], "benchmark_data"])
      String sb_prefix = output_prefix + "." + bd_name

      Array[Array[File?]] bd_snv_ppv_beds = if transpose_site_benchmarking_nested_arrays 
                                            then transpose(site_benchmark_common_snv_ppv_beds[site_bench_di])
                                            else site_benchmark_common_snv_ppv_beds[site_bench_di]
      Array[Array[File?]] bd_indel_ppv_beds = if transpose_site_benchmarking_nested_arrays 
                                              then transpose(site_benchmark_common_indel_ppv_beds[site_bench_di])
                                              else site_benchmark_common_indel_ppv_beds[site_bench_di]
      Array[Array[File?]] bd_sv_ppv_beds = if transpose_site_benchmarking_nested_arrays 
                                           then transpose(site_benchmark_common_sv_ppv_beds[site_bench_di])
                                           else site_benchmark_common_sv_ppv_beds[site_bench_di]
      Array[Array[File?]] bd_snv_sens_beds = if transpose_site_benchmarking_nested_arrays 
                                             then transpose(site_benchmark_common_snv_sens_beds[site_bench_di])
                                             else site_benchmark_common_snv_sens_beds[site_bench_di]
      Array[Array[File?]] bd_indel_sens_beds = if transpose_site_benchmarking_nested_arrays 
                                               then transpose(site_benchmark_common_indel_sens_beds[site_bench_di])
                                               else site_benchmark_common_indel_sens_beds[site_bench_di]
      Array[Array[File?]] bd_sv_sens_beds = if transpose_site_benchmarking_nested_arrays 
                                            then transpose(site_benchmark_common_sv_sens_beds[site_bench_di])
                                            else site_benchmark_common_sv_sens_beds[site_bench_di]
      Array[Array[File?]] bd_ppv_by_af = if transpose_site_benchmarking_nested_arrays 
                                         then transpose(site_benchmark_ppv_by_freqs[site_bench_di])
                                         else site_benchmark_ppv_by_freqs[site_bench_di]
      Array[Array[File?]] bd_sens_by_af = if transpose_site_benchmarking_nested_arrays 
                                          then transpose(site_benchmark_sensitivity_by_freqs[site_bench_di])
                                          else site_benchmark_sensitivity_by_freqs[site_bench_di]

      call PSB.PrepSiteBenchDataToPlot as PrepSiteBench {
        input:
          dataset_name = bd_name,
          dataset_prefix = sb_prefix,
          interval_set_names = select_all(benchmark_interval_names),
          common_snv_ppv_beds = bd_snv_ppv_beds,
          common_indel_ppv_beds = bd_indel_ppv_beds,
          common_sv_ppv_beds = bd_sv_ppv_beds,
          ppv_by_freqs = bd_ppv_by_af,
          sensitivity_by_freqs = bd_sens_by_af,
          bcftools_docker = bcftools_docker,
          g2c_analysis_docker = g2c_analysis_docker
      }

      # Also prep sensitivity benchmarking data
      # This is required to spike in false negative sites into AF benchmarking plots
      call PSB.PrepSiteBenchDataToPlot as PrepSiteBenchInverted {
        input:
          dataset_name = bd_name,
          dataset_prefix = sb_prefix,
          interval_set_names = select_all(benchmark_interval_names),
          common_snv_ppv_beds = bd_snv_sens_beds,
          common_indel_ppv_beds = bd_indel_sens_beds,
          common_sv_ppv_beds = bd_sv_sens_beds,
          bcftools_docker = bcftools_docker,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
  }

  # Preprocess sample benchmarking, if provided
  if (has_sample_benchmarking) {
    scatter ( sample_bench_di in range(n_gb_datasets) ) {

      String gbd_name = select_first([sample_benchmark_dataset_prefixes[sample_bench_di], "benchmark_data"])
      String gb_prefix = output_prefix + "." + gbd_name

      Array[Array[File?]] gbd_ppv_tsvs = if transpose_sample_benchmarking_nested_arrays
                                         then transpose(sample_benchmark_ppv_distribs[sample_bench_di])
                                         else sample_benchmark_ppv_distribs[sample_bench_di]
      Array[Array[File?]] gbd_sens_tsvs = if transpose_sample_benchmarking_nested_arrays
                                          then transpose(sample_benchmark_sensitivity_distribs[sample_bench_di])
                                          else sample_benchmark_sensitivity_distribs[sample_bench_di]

      call PSB.PrepSiteBenchDataToPlot as PrepSampleBench {
        input:
          dataset_name = gbd_name,
          dataset_prefix = gb_prefix,
          interval_set_names = select_all(benchmark_interval_names),
          ppv_by_freqs = gbd_ppv_tsvs,
          sensitivity_by_freqs = gbd_sens_tsvs,
          n_key_columns = 5,
          bcftools_docker = bcftools_docker,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
  }

  # Preprocess twin benchmarking, if provided
  if (has_twin_benchmarking) {
    scatter ( twin_int_i in range(n_bench_intervals) ) {

      Array[File] tb_array_presum = select_all(twin_genotype_benchmark_distribs[twin_int_i])
      String tb_prefix = select_first([select_all(benchmark_interval_names)[twin_int_i], "benchmarking_intervals"])

      call QcTasks.SumCompressedDistribs as SumTwinBenchDistribs {
        input:
          distrib_tsvs = tb_array_presum,
          out_prefix = "twin_bench." + tb_prefix + ".concordance_distribution",
          n_key_columns = 5,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
  }
  Array[File] twin_bench_summed_distribs = select_first([SumTwinBenchDistribs.merged_distrib, []])

  # Preprocess trio benchmarking, if provided
  if (has_trio_benchmarking) {
    scatter ( trio_int_i in range(n_bench_intervals) ) {

      Array[File] mtb_array_presum = select_all(trio_mendelian_violation_distribs[trio_int_i])
      String mtb_prefix = select_first([select_all(benchmark_interval_names)[trio_int_i], "benchmarking_intervals"])

      call QcTasks.SumCompressedDistribs as SumTrioBenchDistribs {
        input:
          distrib_tsvs = mtb_array_presum,
          out_prefix = "trio_bench." + mtb_prefix + ".concordance_distribution",
          n_key_columns = 4,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
  }
  Array[File] trio_bench_summed_distribs = select_first([SumTrioBenchDistribs.merged_distrib, []])


  #################
  ### VISUALIZATION
  #################

  # Generate site-level plots
  call PlotSiteMetrics {
    input:
      size_distrib = size_distrib,
      af_distrib = af_distrib,
      joint_distrib = joint_distrib,
      common_af_cutoff = common_af_cutoff,
      ref_size_distrib = ref_size_distrib,
      ref_af_distrib = ref_af_distrib,
      all_svs_bed = all_svs_bed,
      common_snvs_bed = common_snvs_bed,
      common_indels_bed = common_indels_bed,
      common_svs_bed = common_svs_bed,
      ld_stats_tsv = ld_stats_tsv,
      output_prefix = output_prefix,
      ref_title = ref_cohort_plot_title,
      custom_plotting_constants = custom_plotting_constants,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Generate sample-level plots, including internal benchmarks (twins, trios)
  call PlotSampleMetrics {
    input:
      gt_distrib = gt_distrib,
      ancestry_labels = sample_ancestry_labels,
      phenotype_labels = sample_phenotype_labels,
      sample_subset_list = sample_subset_list,
      twin_concordance_tsvs = twin_bench_summed_distribs,
      trio_concordance_tsvs = trio_bench_summed_distribs,
      eval_interval_names = select_all(benchmark_interval_names),
      common_af_cutoff = common_af_cutoff,
      output_prefix = output_prefix,
      custom_plotting_constants = custom_plotting_constants,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Plot site benchmarking, if provided
  if (has_site_benchmarking) {
    scatter ( site_bench_di in range(n_sb_datasets) ) {

      String plot_bd_prefix = select_all(site_benchmark_dataset_prefixes)[site_bench_di]
      String plot_bd_title = select_all(site_benchmark_dataset_titles)[site_bench_di]

      # Generate plots
      call PlotSiteBenchmarking {
        input:
          ref_dataset_prefix = plot_bd_prefix,
          ref_dataset_title = plot_bd_title,
          eval_interval_names = select_all(benchmark_interval_names),
          inputs_json = select_first([PrepSiteBench.plot_files_json])[site_bench_di],
          inverted_inputs_json = select_first([PrepSiteBenchInverted.plot_files_json])[site_bench_di],
          output_prefix = output_prefix,
          common_af_cutoff = common_af_cutoff,
          custom_plotting_constants = custom_plotting_constants,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
  }

  # Plot sample benchmarking, if provided
  if (has_sample_benchmarking) {
    scatter ( sample_bench_di in range(n_gb_datasets) ) {

      String plot_gbd_prefix = select_all(sample_benchmark_dataset_prefixes)[sample_bench_di]
      String plot_gbd_title = select_all(sample_benchmark_dataset_titles)[sample_bench_di]

      # Generate plots
      call PlotSampleBenchmarking {
        input:
          ref_dataset_prefix = plot_gbd_prefix,
          ref_dataset_title = plot_gbd_title,
          eval_interval_names = select_all(benchmark_interval_names),
          inputs_json = select_first([PrepSampleBench.plot_files_json])[sample_bench_di],
          sample_subset_list = sample_subset_list,
          output_prefix = output_prefix,
          common_af_cutoff = common_af_cutoff,
          custom_plotting_constants = custom_plotting_constants,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
  }


  ##################
  ### OUTPUT CLEANUP
  ##################

  # Package all outputs into two tarballs; one for plots and one for stats
  Array[Array[File]] all_out_tarballs_preflat = [[PlotSiteMetrics.site_metric_plots_tarball],
                                                 [PlotSampleMetrics.sample_metric_plots_tarball],
                                                 select_first([PlotSiteBenchmarking.site_benchmarking_plots_tarball, []]),
                                                 select_first([PlotSampleBenchmarking.sample_benchmarking_plots_tarball, []])]
  Array[File] all_out_tarballs = flatten(select_all(all_out_tarballs_preflat))
  call PackageOutputs {
    input:
      tarballs = all_out_tarballs,
      previous_stats = previous_stats,
      site_benchmark_prefix = select_all(site_benchmark_dataset_prefixes),
      site_benchmark_title = select_all(site_benchmark_dataset_titles),
      sample_benchmark_prefixes = select_all(sample_benchmark_dataset_prefixes),
      sample_benchmark_titles = select_all(sample_benchmark_dataset_titles),
      custom_targets = custom_qc_target_metrics,
      out_prefix = output_prefix,
      custom_plotting_constants = custom_plotting_constants,
      g2c_analysis_docker = g2c_analysis_docker
  }

  output {
    File plots_tarball = PackageOutputs.plots_tarball
    File stats_tarball = PackageOutputs.stats_tarball
    File all_stats_tsv = PackageOutputs.all_stats_tsv
  }
}


task PackageOutputs {
  input {
    Array[File] tarballs
    File? previous_stats
    Array[String?] site_benchmark_prefix
    Array[String?] site_benchmark_title
    Array[String?] sample_benchmark_prefixes
    Array[String?] sample_benchmark_titles
    File? custom_targets
    String out_prefix
    File? custom_plotting_constants
    String g2c_analysis_docker
  }

  String previous_opt = if defined(previous_stats) then "--previous-stats ~{basename(select_first([previous_stats, 'not_real.txt']))}" else ""
  String targets_opt = if defined(custom_targets) then "--custom-targets ~{basename(select_first([custom_targets, 'not_real.txt']))}" else ""
  String constants_opt = if defined(custom_plotting_constants) then "--custom-constants ~{basename(select_first([custom_plotting_constants, 'not_real.txt']))}"  else ""

  Int disk_gb = ceil(10 * size(tarballs, "GB")) + 10

  command <<<
    set -eu -o pipefail

    # Make output collection directories
    for subset in plots stats; do
      mkdir ~{out_prefix}.$subset
    done
    mkdir ~{out_prefix}.plots/~{out_prefix}.qc_summary

    # Unpack all plots
    while read tb; do
      tar -xzvf $tb -C ~{out_prefix}.plots/
    done < ~{write_lines(tarballs)}

    # Move all .tsvs and .txt files to stats
    for suffix in tsv tsv.gz txt txt.gz; do
      find ~{out_prefix}.plots -name "*.$suffix" || true
    done | xargs -I {} mv {} ~{out_prefix}.stats/ || true

    # Link previous stats and custom targets to pwd if provided
    if ~{defined(previous_stats)}; then
      ln -s ~{default="" previous_stats} ./
    fi
    if ~{defined(custom_targets)}; then
      ln -s ~{default="" custom_targets} ./
    fi
    if ~{defined(custom_plotting_constants)}; then
      ln -s ~{default="" custom_plotting_constants} ./
    fi

    # Collapse all summary stats and generate summary barplots
    echo -e "#analysis\tmeasure\tvalue\tn" > ss.all.tsv
    find ~{out_prefix}.stats/ -name "*.summary_metrics.tsv" \
    | xargs -I {} cat {} | grep -ve '^#' | grep -ve '^analysis' \
    | sort -Vk1,1 -k2,2V -k3,3n -k4,4n >> \
    ~{out_prefix}.all_qc_summary_metrics.tsv
    cp \
      ~{out_prefix}.all_qc_summary_metrics.tsv \
      ~{out_prefix}.stats/~{out_prefix}.all_qc_summary_metrics.tsv
    cmd="Rscript /opt/pancan_germline_wgs/scripts/qc/vcf_qc/plot_overall_qc_summary.R"
    cmd="$cmd --stats ~{out_prefix}.stats/~{out_prefix}.all_qc_summary_metrics.tsv"
    cmd="$cmd ~{previous_opt} ~{targets_opt}"
    while read sbp; do
      cmd="$cmd --site-ref-prefix \"$sbp\""
    done < ~{write_lines(select_all(site_benchmark_prefix))}
    while read sbt; do
      cmd="$cmd --site-ref-title \"$sbt\""
    done < ~{write_lines(select_all(site_benchmark_title))}
    while read gbp; do
      cmd="$cmd --sample-benchmarking-prefix \"$gbp\""
    done < ~{write_lines(select_all(sample_benchmark_prefixes))}
    while read gbt; do
      cmd="$cmd --sample-benchmarking-title \"$gbt\""
    done < ~{write_lines(select_all(sample_benchmark_titles))}
    cmd="$cmd ~{constants_opt} --out-prefix \"~{out_prefix}.plots/~{out_prefix}.qc_summary/~{out_prefix}\""
    echo -e "Now generating summary plots as follows:\n\n$cmd"
    eval "$cmd"

    # Compress outputs
    for subset in plots stats; do
      tar -czvf ~{out_prefix}.$subset.tar.gz ~{out_prefix}.$subset
    done
  >>>

  output {
    File plots_tarball = "~{out_prefix}.plots.tar.gz"
    File stats_tarball = "~{out_prefix}.stats.tar.gz"
    File all_stats_tsv = "~{out_prefix}.all_qc_summary_metrics.tsv"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 2
    max_retries: 1
  }
}


task PlotSampleBenchmarking {
  input {
    String ref_dataset_prefix
    String ref_dataset_title
    Array[String] eval_interval_names
    File inputs_json
    File? sample_subset_list
    String output_prefix
    Float common_af_cutoff
    File? custom_plotting_constants

    Float mem_gb = 7.5
    Int n_cpu = 4
    Int disk_gb = 50

    String g2c_analysis_docker
  }

  String subset_opt = if defined(sample_subset_list) then "--subset-samples ~{basename(select_first([sample_subset_list, 'not_real.txt']))}"  else ""
  String constants_opt = if defined(custom_plotting_constants) then "--custom-constants ~{basename(select_first([custom_plotting_constants, 'not_real.txt']))}"  else ""

  # Note that this outdir string is used as both a directory name and a file prefix
  String outdir = sub(output_prefix + "." + ref_dataset_prefix + "." + "sample_benchmarking", "[ ]+", "_")

  Int n_sets = length(eval_interval_names)

  command <<<
    set -euo pipefail

    mkdir ~{outdir}

    # Parse inputs .json generated by PSB
    python3 - "~{inputs_json}" <<CODE
import json
import sys

infile = sys.argv[1]

with open(infile, "r") as f:
  data = json.load(f)

def write_array(name, val):
  val = [v for v in val if v is not None and v != ""]
  if len(val) == 0:
    return
  with open(f"{name}.txt", "w") as out:
      for item in (val or []):
          if item is not None:
              out.write(item + "\n")

write_array("ppv_by_af", data.get("ppv_by_af"))
write_array("sens_by_af", data.get("sens_by_af"))
CODE

    # Relocate input files
    if ~{defined(sample_subset_list)}; then
      ln -s ~{default="" sample_subset_list} ./
    fi
    if ~{defined(custom_plotting_constants)}; then
      ln -s ~{default="" custom_plotting_constants} ./
    fi

    # Plot sample summary metrics like PPV and sensitivity
    cmd="/opt/pancan_germline_wgs/scripts/qc/vcf_qc/plot_external_sample_benchmarking.R"
    if [ -s ppv_by_af.txt ]; then
      cat ppv_by_af.txt | gsutil -m cp -I ./
      while read tsv; do
        cmd="$cmd --ppv-tsv $tsv"
      done < <( cat ppv_by_af.txt | xargs -I {} basename {} )
    fi
    if [ -s sens_by_af.txt ]; then
      cat sens_by_af.txt | gsutil -m cp -I ./
      while read tsv; do
        cmd="$cmd --sens-tsv $tsv"
      done < <( cat sens_by_af.txt | xargs -I {} basename {} )
    fi
    while read sname; do
      cmd="$cmd --set-name $sname"
    done < ~{write_lines(eval_interval_names)}
    cmd="$cmd --ref-title \"~{ref_dataset_title}\" --common-af ~{common_af_cutoff} ~{subset_opt}"
    cmd="$cmd ~{constants_opt} --out-prefix ~{outdir}/~{outdir}"
    echo -e "Now performing sample benchmarking metric visualization as follows:\n$cmd"
    eval "$cmd"

    # Compress outputs
    tar -czvf "~{outdir}.tar.gz" ~{outdir}/
  >>>

  output {
    File sample_benchmarking_plots_tarball = "~{outdir}.tar.gz"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 2
    max_retries: 1
  }
}


task PlotSampleMetrics {
  input {
    File gt_distrib
    File? ancestry_labels
    File? phenotype_labels
    File? sample_subset_list

    Array[File] twin_concordance_tsvs = []
    Array[File] trio_concordance_tsvs = []
    Array[String] eval_interval_names = []
    
    Float common_af_cutoff
    String output_prefix
    File? custom_plotting_constants

    Float mem_gb = 7.5
    Int n_cpu = 4
    Int? disk_gb

    String g2c_analysis_docker
  }

  String pop_opt = if defined(ancestry_labels) then "--ancestry-labels ~{basename(select_first([ancestry_labels, 'not_real.txt']))}" else ""
  String pheno_opt = if defined(phenotype_labels) then "--phenotype-labels ~{basename(select_first([phenotype_labels, 'not_real.txt']))}"  else ""
  String subset_opt = if defined(sample_subset_list) then "--subset-samples ~{basename(select_first([sample_subset_list, 'not_real.txt']))}"  else ""
  String constants_opt = if defined(custom_plotting_constants) then "--custom-constants ~{basename(select_first([custom_plotting_constants, 'not_real.txt']))}"  else ""

  Boolean do_twins = length(eval_interval_names) > 0 && length(twin_concordance_tsvs) > 0
  Boolean do_trios = length(eval_interval_names) > 0 && length(trio_concordance_tsvs) > 0

  Int default_disk_gb = ceil(2 * size(flatten([[gt_distrib], twin_concordance_tsvs, trio_concordance_tsvs]), "GB")) + 20

  command <<<
    set -eu -o pipefail

    mkdir ~{output_prefix}.sample_metrics

    # Relocate optional inputs where necessary
    if ~{defined(ancestry_labels)}; then
      ln -s ~{default="" ancestry_labels} ./
    fi
    if ~{defined(phenotype_labels)}; then
      ln -s ~{default="" phenotype_labels} ./
    fi
    if ~{defined(sample_subset_list)}; then
      ln -s ~{default="" sample_subset_list} ./
    fi
    if ~{defined(custom_plotting_constants)}; then
      ln -s ~{default="" custom_plotting_constants} ./
    fi

    # Plot variation per genome
    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/plot_variants_per_sample.R \
      --genotype-dist-tsv ~{gt_distrib} \
      ~{pop_opt} \
      ~{pheno_opt} \
      ~{subset_opt} \
      --out-prefix ~{output_prefix}.sample_metrics/~{output_prefix}

    # Plot twin benchmarking, if optioned
    if ~{do_twins}; then
      twin_cmd="/opt/pancan_germline_wgs/scripts/qc/vcf_qc/plot_twin_benchmarking.R"
      while read tsv; do
        twin_cmd="$twin_cmd --bench-tsv $tsv"
      done < ~{write_lines(twin_concordance_tsvs)}
      while read sname; do
        twin_cmd="$twin_cmd --set-name \"$sname\""
      done < ~{write_lines(eval_interval_names)}
      twin_cmd="$twin_cmd --common-af ~{common_af_cutoff} ~{subset_opt}"
      twin_cmd="$twin_cmd ~{constants_opt} --out-prefix ~{output_prefix}.sample_metrics/~{output_prefix}"
      echo -e "Now performing twin benchmarking as follows:\n$twin_cmd"
      eval "$twin_cmd"
    fi

    # Plot trio benchmarking, if optioned
    if ~{do_trios}; then
      trio_cmd="/opt/pancan_germline_wgs/scripts/qc/vcf_qc/plot_trio_benchmarking.R"
      while read tsv; do
        trio_cmd="$trio_cmd --bench-tsv $tsv"
      done < ~{write_lines(trio_concordance_tsvs)}
      while read sname; do
        trio_cmd="$trio_cmd --set-name \"$sname\""
      done < ~{write_lines(eval_interval_names)}
      trio_cmd="$trio_cmd --common-af ~{common_af_cutoff}"
      trio_cmd="$trio_cmd ~{constants_opt} --out-prefix ~{output_prefix}.sample_metrics/~{output_prefix}"
      echo -e "Now performing trio benchmarking as follows:\n$trio_cmd"
      eval "$trio_cmd"
    fi

    # Compress outputs
    tar -czvf ~{output_prefix}.sample_metrics.tar.gz ~{output_prefix}.sample_metrics/
  >>>

  output {
    File sample_metric_plots_tarball = "~{output_prefix}.sample_metrics.tar.gz"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    preemptible: 2
    max_retries: 1
  }
}


task PlotSiteBenchmarking {
  input {
    String ref_dataset_prefix
    String ref_dataset_title
    Array[String] eval_interval_names
    Boolean bed_arrays_include_union = true
    File inputs_json
    File inverted_inputs_json
    String output_prefix
    Float common_af_cutoff
    File? custom_plotting_constants

    Float mem_gb = 7.5
    Int n_cpu = 4
    Int disk_gb = 50

    String g2c_analysis_docker
  }

  String constants_opt = if defined(custom_plotting_constants) then "--custom-constants ~{basename(select_first([custom_plotting_constants, 'not_real.txt']))}"  else ""

  # Note that this outdir string is used as both a directory name and a file prefix
  String outdir = sub(output_prefix + "." + ref_dataset_prefix + "." + "site_benchmarking", "[ ]+", "_")

  Array[String] eval_interval_names_plus_union = if bed_arrays_include_union 
                                                 then flatten([eval_interval_names, ["All"]]) 
                                                 else eval_interval_names
  Int n_sets = length(eval_interval_names_plus_union)

  command <<<
    set -euo pipefail

    mkdir ~{outdir}

    # Parse inputs .json generated by PSB
    python3 - "~{inputs_json}" <<CODE
import json
import sys

infile = sys.argv[1]
with open(infile, "r") as f:
  data = json.load(f)

def write_array(name, val):
  val = [v for v in val if v is not None and v != ""]
  if len(val) == 0:
    return
  with open(f"{name}.txt", "w") as out:
      for item in (val or []):
          if item is not None:
              out.write(item + "\n")

write_array("snv_ppv_beds", data.get("snv_ppv_beds", []) + [data.get("snv_ppv_union")])
write_array("indel_ppv_beds", data.get("indel_ppv_beds", []) + [data.get("indel_ppv_union")])
write_array("sv_ppv_beds", data.get("sv_ppv_beds", []) + [data.get("sv_ppv_union")])

write_array("ppv_by_af", data.get("ppv_by_af"))
write_array("sens_by_af", data.get("sens_by_af"))
CODE

    # Parse inputs .json for "inverted" files (to spike false negatives into plots)
    python3 - "~{inverted_inputs_json}" <<CODE
import json
import sys

infile = sys.argv[1]
with open(infile, "r") as f:
  data = json.load(f)

def write_array(name, val):
  val = [v for v in val if v is not None and v != ""]
  if len(val) == 0:
    return
  with open(f"{name}.txt", "w") as out:
      for item in (val or []):
          if item is not None:
              out.write(item + "\n")

write_array("snv_sens_beds", data.get("snv_ppv_beds", []) + [data.get("snv_ppv_union")])
write_array("indel_sens_beds", data.get("indel_ppv_beds", []) + [data.get("indel_ppv_union")])
write_array("sv_sens_beds", data.get("sv_ppv_beds", []) + [data.get("sv_ppv_union")])
CODE

    # Localize SNV beds
    if [ -s snv_ppv_beds.txt ]; then
      cat snv_ppv_beds.txt | gsutil -m cp -I ./
      cat snv_ppv_beds.txt | xargs -I {} basename {} > snv_beds.list
      # Add false negatives to BED for plotting
      if [ -s snv_sens_beds.txt ]; then
        while read main_bed sens_uri; do
          gsutil -m cat $sens_uri \
          | gunzip -c \
          | awk -v FS="\t" -v OFS="\t" \
            '{ if ($9=="NA") print $1, $2, $3, "FN_"$4, $5, $6, $7, $10, $9, $8, $11 }' \
          | cat <( zcat $main_bed ) - \
          | gzip -c \
          > tmp.bed.gz
          mv tmp.bed.gz $main_bed
        done < <( paste snv_beds.list snv_sens_beds.txt )
      fi
    else
      seq 1 ~{n_sets} | awk '{ print "." }' > snv_beds.list
    fi

    # Localize indel beds
    if [ -s indel_ppv_beds.txt ]; then
      cat indel_ppv_beds.txt | gsutil -m cp -I ./
      cat indel_ppv_beds.txt | xargs -I {} basename {} > indel_beds.list
      # Add false negatives to BED for plotting
      if [ -s indel_sens_beds.txt ]; then
        while read main_bed sens_uri; do
          gsutil -m cat $sens_uri \
          | gunzip -c \
          | awk -v FS="\t" -v OFS="\t" \
            '{ if ($9=="NA") print $1, $2, $3, "FN_"$4, $5, $6, $7, $10, $9, $8, $11 }' \
          | cat <( zcat $main_bed ) - \
          | gzip -c \
          > tmp.bed.gz
          mv tmp.bed.gz $main_bed
        done < <( paste indel_beds.list indel_sens_beds.txt )
      fi
    else
      seq 1 ~{n_sets} | awk '{ print "." }' > indel_beds.list
    fi

    # Localize SV beds
    if [ -s sv_ppv_beds.txt ]; then
      cat sv_ppv_beds.txt | gsutil -m cp -I ./
      cat sv_ppv_beds.txt | xargs -I {} basename {} > sv_beds.list
      # Add false negatives to BED for plotting
      if [ -s sv_sens_beds.txt ]; then
        while read main_bed sens_uri; do
          gsutil -m cat $sens_uri \
          | gunzip -c \
          | awk -v FS="\t" -v OFS="\t" \
            '{ if ($9=="NA") print $1, $2, $3, "FN_"$4, $5, $6, $7, $10, $9, $8, $11 }' \
          | cat <( zcat $main_bed ) - \
          | gzip -c \
          > tmp.bed.gz
          mv tmp.bed.gz $main_bed
        done < <( paste sv_beds.list sv_sens_beds.txt )
      fi
    else
      seq 1 ~{n_sets} | awk '{ print "." }' > sv_beds.list
    fi

    # Relocate custom constants if provided
    if ~{defined(custom_plotting_constants)}; then
      ln -s ~{default="" custom_plotting_constants} ./
    fi

    # Make input .tsv for pointwise plotting
    paste \
      ~{write_lines(eval_interval_names_plus_union)} \
      snv_beds.list \
      indel_beds.list \
      sv_beds.list \
    > pointwise.in.tsv
    echo "Pointwise benchmarking inputs parsed as follows:"
    cat pointwise.in.tsv

    # Perform pointwise plotting once for each eval interval set
    while read set_name snv_bed indel_bed sv_bed; do
      cmd="/opt/pancan_germline_wgs/scripts/qc/vcf_qc/plot_pointwise_site_benchmarking.R"
      if [ $snv_bed != "." ]; then
        cmd="$cmd --snvs $snv_bed"
      fi
      if [ $indel_bed != "." ]; then
        cmd="$cmd --indels $indel_bed"
      fi
      if [ $sv_bed != "." ]; then
        cmd="$cmd --svs $sv_bed"
      fi
      cmd="$cmd --common-af ~{common_af_cutoff} --ref-title \"~{ref_dataset_title}\""
      cmd="$cmd --combine ~{constants_opt} --out-prefix ~{outdir}/~{outdir}.$set_name --set-name $set_name"
      echo -e "Now performing pointwise site benchmarking as follows:\n$cmd"
      eval "$cmd"
    done < pointwise.in.tsv

    # Plot site summary metrics like PPV and sensitivity
    cmd="/opt/pancan_germline_wgs/scripts/qc/vcf_qc/plot_site_benchmarking_metrics.R"
    if [ -s ppv_by_af.txt ]; then
      cat ppv_by_af.txt | gsutil -m cp -I ./
      while read tsv; do
        cmd="$cmd --ppv-by-af $tsv"
      done < <( cat ppv_by_af.txt | xargs -I {} basename {} )
    fi
    if [ -s sens_by_af.txt ]; then
      cat sens_by_af.txt | gsutil -m cp -I ./
      while read tsv; do
        cmd="$cmd --sens-by-af $tsv"
      done < <( cat sens_by_af.txt | xargs -I {} basename {} )
    fi
    while read sname; do
      cmd="$cmd --set-name $sname"
    done < ~{write_lines(eval_interval_names)}
    cmd="$cmd --ref-title \"~{ref_dataset_title}\" --common-af ~{common_af_cutoff}"
    cmd="$cmd ~{constants_opt} --out-prefix ~{outdir}/~{outdir}"
    echo -e "Now performing site benchmarking metric visualization as follows:\n$cmd"
    eval "$cmd"

    # Compress outputs
    tar -czvf "~{outdir}.tar.gz" ~{outdir}/
  >>>

  output {
    File site_benchmarking_plots_tarball = "~{outdir}.tar.gz"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 2
    max_retries: 1
  }
}


task PlotSiteMetrics {
  input {
    File size_distrib
    File af_distrib
    File joint_distrib
    
    Float common_af_cutoff
    
    File? ref_size_distrib
    File? ref_af_distrib
    
    File? all_svs_bed
    File? common_snvs_bed
    File? common_indels_bed
    File? common_svs_bed

    File? ld_stats_tsv

    String output_prefix
    String? ref_title
    File? custom_plotting_constants

    Float mem_gb = 7.5
    Int n_cpu = 4
    Int? disk_gb

    String g2c_analysis_docker
  }

  Array[File?] loc_inputs = [size_distrib, af_distrib, joint_distrib, all_svs_bed,
                             common_snvs_bed, common_indels_bed, common_svs_bed, 
                             ld_stats_tsv]
  Int default_disk_gb = ceil(2 * size(select_all(loc_inputs), "GB")) + 20

  String constants_opt = if defined(custom_plotting_constants) then "--custom-constants ~{basename(select_first([custom_plotting_constants, 'not_real.txt']))}"  else ""

  Boolean has_ref_size = defined(ref_size_distrib)
  String ref_size_bname = if has_ref_size then basename(select_first([ref_size_distrib])) else ""
  String ref_size_cmd = if has_ref_size then "--ref-size-distrib ~{ref_size_bname}" else ""

  Boolean has_ref_af = defined(ref_af_distrib)
  String ref_af_bname = if has_ref_af then basename(select_first([ref_af_distrib])) else ""
  String ref_af_cmd = if has_ref_af then "--ref-af-distrib ~{ref_af_bname}" else ""

  String ref_title_cmd = if defined(ref_title) then "--ref-title \"~{ref_title}\"" else ""

  Boolean has_all_svs = defined(all_svs_bed)
  String all_sv_bname = if has_all_svs then basename(select_first([all_svs_bed])) else ""
  String summary_sv_cmd = if has_all_svs then "--sv-sites ~{all_sv_bname}" else ""

  Boolean has_common_snvs = defined(common_snvs_bed)
  String common_snv_bname = if has_common_snvs then basename(select_first([common_snvs_bed])) else ""
  String pw_snv_cmd = if has_common_snvs then "--snvs ~{common_snv_bname}" else ""
  
  Boolean has_common_indels = defined(common_indels_bed)
  String common_indel_bname = if has_common_indels then basename(select_first([common_indels_bed])) else ""
  String pw_indel_cmd = if has_common_indels then "--indels ~{common_indel_bname}" else ""

  Boolean has_common_svs = defined(common_svs_bed)
  String common_sv_bname = if has_common_svs then basename(select_first([common_svs_bed])) else ""
  String pw_sv_cmd = if has_common_svs then "--svs ~{common_sv_bname}" else ""

  Boolean has_ld = defined(ld_stats_tsv)
  String ld_stat_bname = if has_ld then basename(select_first([ld_stats_tsv])) else ""
  String ld_cmd = if has_ld then "--ld-stats ~{ld_stat_bname}" else ""

  Boolean has_common_variants = has_common_snvs || has_common_indels || has_common_svs

  command <<<
    set -eu -o pipefail

    mkdir ~{output_prefix}.site_metrics

    # Symlink optional inputs to working directory
    if ~{has_all_svs}; then
      ln -s ~{default="" all_svs_bed} ~{all_sv_bname}
    fi
    if ~{has_ref_size}; then
      ln -s ~{default="" ref_size_distrib} ~{ref_size_bname}
    fi
    if ~{has_ref_af}; then
      ln -s ~{default="" ref_af_distrib} ~{ref_af_bname}
    fi
    if ~{defined(custom_plotting_constants)}; then
      ln -s ~{default="" custom_plotting_constants} ./
    fi

    # Plot site summary metrics
    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/plot_site_summary_metrics.R \
      --size-distrib ~{size_distrib} \
      --af-distrib ~{af_distrib} \
      --joint-distrib ~{joint_distrib} \
      ~{ref_size_cmd} \
      ~{ref_af_cmd} \
      ~{ref_title_cmd} \
      ~{summary_sv_cmd} \
      --common-af ~{common_af_cutoff} \
      ~{constants_opt} \
      --out-prefix ~{output_prefix}.site_metrics/~{output_prefix}

    # Symlink common variant BEDs to working directory
    if ~{has_common_snvs}; then
      ln -s ~{default="" common_snvs_bed} ~{common_snv_bname}
    fi
    if ~{has_common_indels}; then
      ln -s ~{default="" common_indels_bed} ~{common_indel_bname}
    fi
    if ~{has_common_svs}; then
      ln -s ~{default="" common_svs_bed} ~{common_sv_bname}
    fi

    # Plot site-level metrics for common variants
    if ~{has_common_variants}; then

      # Symlink LD stats to working directory, if provided
      if ~{has_ld}; then
        ln -s ~{default="" ld_stats_tsv} ~{ld_stat_bname}
      fi

      # Plot pointwise QC
      /opt/pancan_germline_wgs/scripts/qc/vcf_qc/plot_site_pointwise_metrics.R \
        ~{pw_snv_cmd} \
        ~{pw_indel_cmd} \
        ~{pw_sv_cmd} \
        --combine \
        --common-af ~{common_af_cutoff} \
        ~{ld_cmd} \
        ~{constants_opt} \
        --out-prefix ~{output_prefix}.site_metrics/~{output_prefix}

      # Append summary stats to previous stats file
      fgrep -v "#" \
        ~{output_prefix}.site_metrics/~{output_prefix}.pointwise_summary_metrics.tsv \
      >> ~{output_prefix}.site_metrics/~{output_prefix}.summary_metrics.tsv
      rm ~{output_prefix}.site_metrics/~{output_prefix}.pointwise_summary_metrics.tsv
    fi

    # Compress outputs
    tar -czvf ~{output_prefix}.site_metrics.tar.gz ~{output_prefix}.site_metrics/
  >>>

  output {
    File site_metric_plots_tarball = "~{output_prefix}.site_metrics.tar.gz"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    preemptible: 2
    max_retries: 1
  }
}

