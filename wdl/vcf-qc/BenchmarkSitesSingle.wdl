# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Perform two-sided site-level benchmarking of a "source" callset versus a reference "target" callset
# Note: this is a subroutine called within BenchmarkSites.wdl, and is designed to be
# invoked by this parent workflow

# Expects that both callsets have already been processed by CollectVcfQcMetrics.wdl


version 1.0


import "QcTasks.wdl" as QcTasks


workflow BenchmarkSitesSingle {
  input {
    File eval_interval_bed
    File genome_file

    File? source_snv_bed
    File? source_snv_bed_idx
    File? target_snv_bed
    File? target_snv_bed_idx

    File? source_indel_bed
    File? source_indel_bed_idx
    File? target_indel_bed
    File? target_indel_bed_idx

    File? source_sv_bed
    File? source_sv_bed_idx
    File? target_sv_bed
    File? target_sv_bed_idx

    String eval_prefix
    String source_prefix
    String target_prefix

    Float common_af_cutoff

    Int total_shards

    Boolean make_full_id_maps = false

    Float compare_sites_mem_gb = 3.75
    Float snv_mem_scalar = 2.0
    Float indel_mem_scalar = 1.0

    Float shard_intervals_mem_gb = 3.75

    String bcftools_docker
    String g2c_analysis_docker
  }

  String ppv_prefix = eval_prefix + "." + source_prefix + "_vs_" + target_prefix
  String sens_prefix = eval_prefix + "." + target_prefix + "_vs_" + source_prefix

  Boolean do_snv = if defined(source_snv_bed) && defined(source_snv_bed_idx)
                   && defined(target_snv_bed) && defined(target_snv_bed_idx)
                   then true else false
  Boolean do_indel = if defined(source_indel_bed) && defined(source_indel_bed_idx)
                     && defined(target_indel_bed) && defined(target_indel_bed_idx)
                     then true else false
  Boolean do_sv = if defined(source_sv_bed) && defined(source_sv_bed_idx)
                  && defined(target_sv_bed) && defined(target_sv_bed_idx) 
                  then true else false
  Int n_var_types = (if do_snv then 1 else 0) + 
                    (if do_indel then 1 else 0) + 
                    (if do_sv then 1 else 0)

  # Shard evaluation intervals, balancing on genomic bp
  Int target_shards = floor(total_shards / (2 * n_var_types))
  Int shard_ints_n_cpu = ceil(shard_intervals_mem_gb / 2)
  call QcTasks.ShardIntervals {
    input:
      intervals_bed = eval_interval_bed,
      n_shards = target_shards,
      prefix = eval_prefix,
      mem_gb = shard_intervals_mem_gb,
      n_cpu = shard_ints_n_cpu,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Process each eval interval in parallel
  scatter ( eval_shard in ShardIntervals.interval_shards ) {

    String shard_prefix = basename(eval_shard, ".bed.gz")

    # Compare SNVs, if provided
    if ( do_snv ) {
      call QcTasks.PrepSites as PrepSourceSnvs {
        input:
          beds = select_all([source_snv_bed]),
          bed_idxs = select_all([source_snv_bed_idx]),
          eval_interval_bed = eval_shard,
          min_size = 0,
          max_size = 0,
          prefix = source_prefix + ".snv." + shard_prefix,
          g2c_analysis_docker = g2c_analysis_docker
      }

      call QcTasks.PrepSites as PrepTargetSnvs {
        input:
          beds = select_all([target_snv_bed]),
          bed_idxs = select_all([target_snv_bed_idx]),
          eval_interval_bed = eval_shard,
          min_size = 0,
          max_size = 0,
          prefix = target_prefix + ".snv." + shard_prefix,
          g2c_analysis_docker = g2c_analysis_docker
      }

      call CompareSites as CompareSnvPpv {
        input:
          query_bed = PrepSourceSnvs.query_bed,
          ref_bed = PrepTargetSnvs.ref_bed,
          genome_file = genome_file,
          mode = "exact",
          common_af = common_af_cutoff,
          prefix = ppv_prefix + ".snv." + shard_prefix,
          mem_gb = snv_mem_scalar * compare_sites_mem_gb,
          g2c_analysis_docker = g2c_analysis_docker
      }

      call CompareSites as CompareSnvSens {
        input:
          query_bed = PrepTargetSnvs.query_bed,
          ref_bed = PrepSourceSnvs.ref_bed,
          genome_file = genome_file,
          mode = "exact",
          common_af = common_af_cutoff,
          degenerate = true,
          prefix = sens_prefix + ".snv." + shard_prefix,
          mem_gb = snv_mem_scalar * compare_sites_mem_gb,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }

    # Compare indels, if provided
    if ( do_indel ) {
      call QcTasks.PrepSites as PrepSourceIndels {
        input:
          beds = select_all([source_indel_bed, source_sv_bed]),
          bed_idxs = select_all([source_indel_bed_idx, source_sv_bed_idx]),
          eval_interval_bed = eval_shard,
          min_size = 1,
          max_size = 49,
          prefix = source_prefix + ".indel." + shard_prefix,
          g2c_analysis_docker = g2c_analysis_docker
      }

      call QcTasks.PrepSites as PrepTargetIndels {
        input:
          beds = select_all([target_indel_bed, target_sv_bed]),
          bed_idxs = select_all([target_indel_bed_idx, target_sv_bed_idx]),
          eval_interval_bed = eval_shard,
          min_size = 1,
          max_size = 49,
          prefix = target_prefix + ".indel." + shard_prefix,
          g2c_analysis_docker = g2c_analysis_docker
      }

      call CompareSites as CompareIndelPpv {
        input:
          query_bed = PrepSourceIndels.query_bed,
          ref_bed = PrepTargetIndels.ref_bed,
          genome_file = genome_file,
          mode = "both",
          common_af = common_af_cutoff,
          prefix = ppv_prefix + ".indel." + shard_prefix,
          mem_gb = indel_mem_scalar * compare_sites_mem_gb,
          g2c_analysis_docker = g2c_analysis_docker
      }

      call CompareSites as CompareIndelSens {
        input:
          query_bed = PrepTargetIndels.query_bed,
          ref_bed = PrepSourceIndels.ref_bed,
          genome_file = genome_file,
          mode = "both",
          common_af = common_af_cutoff,
          degenerate = true,
          prefix = sens_prefix + ".indel." + shard_prefix,
          mem_gb = indel_mem_scalar * compare_sites_mem_gb,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }

    # Compare SVs, if provided
    if ( do_sv ) {
      call QcTasks.PrepSites as PrepSourceSvs {
        input:
          beds = select_all([source_sv_bed, source_indel_bed]),
          bed_idxs = select_all([source_sv_bed_idx, source_indel_bed_idx]),
          eval_interval_bed = eval_shard,
          min_size = 50,
          max_size = 1000000000,
          prefix = source_prefix + ".sv." + shard_prefix,
          g2c_analysis_docker = g2c_analysis_docker
      }

      call QcTasks.PrepSites as PrepTargetSvs {
        input:
          beds = select_all([target_sv_bed, target_indel_bed]),
          bed_idxs = select_all([target_sv_bed_idx, target_indel_bed_idx]),
          eval_interval_bed = eval_shard,
          min_size = 50,
          max_size = 1000000000,
          prefix = target_prefix + ".sv." + shard_prefix,
          g2c_analysis_docker = g2c_analysis_docker
      }

      call CompareSites as CompareSvPpv {
        input:
          query_bed = PrepSourceSvs.query_bed,
          ref_bed = PrepTargetSvs.ref_bed,
          genome_file = genome_file,
          mode = "both",
          overlap_pad = 2,
          common_af = common_af_cutoff,
          prefix = ppv_prefix + ".sv." + shard_prefix,
          mem_gb = compare_sites_mem_gb,
          g2c_analysis_docker = g2c_analysis_docker
      }

      call CompareSites as CompareSvSens {
        input:
          query_bed = PrepTargetSvs.query_bed,
          ref_bed = PrepSourceSvs.ref_bed,
          genome_file = genome_file,
          mode = "both",
          overlap_pad = 2,
          common_af = common_af_cutoff,
          degenerate = true,
          prefix = sens_prefix + ".sv." + shard_prefix,
          mem_gb = compare_sites_mem_gb,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
  }

  # Collapse SNV common site comparisons
  Array[File] common_snv_ppv_beds = select_all(CompareSnvPpv.common_sites_bed)
  if ( length(common_snv_ppv_beds) > 0 ) {
    call QcTasks.ConcatTextFiles as CollapseCommonSnvPpvSites {
      input:
        shards = common_snv_ppv_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = ppv_prefix + ".common_sites.snvs.bed.gz",
        mem_gb = 3.75,
        n_cpu = 2,
        docker = bcftools_docker
      }
  }
  Array[File] common_snv_sens_beds = select_all(CompareSnvSens.common_sites_bed)
  if ( length(common_snv_sens_beds) > 0 ) {
    call QcTasks.ConcatTextFiles as CollapseCommonSnvSensSites {
      input:
        shards = common_snv_sens_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = sens_prefix + ".common_sites.snvs.bed.gz",
        mem_gb = 3.75,
        n_cpu = 2,
        docker = bcftools_docker
      }
  }

  # Collapse indel common site comparisons
  Array[File] common_indel_ppv_beds = select_all(CompareIndelPpv.common_sites_bed)
  if ( length(common_indel_ppv_beds) > 0 ) {
    call QcTasks.ConcatTextFiles as CollapseCommonIndelPpvSites {
      input:
        shards = common_indel_ppv_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = ppv_prefix + ".common_sites.indels.bed.gz",
        mem_gb = 3.75,
        n_cpu = 2,
        docker = bcftools_docker
      }
  }
  Array[File] common_indel_sens_beds = select_all(CompareIndelSens.common_sites_bed)
  if ( length(common_indel_sens_beds) > 0 ) {
    call QcTasks.ConcatTextFiles as CollapseCommonIndelSensSites {
      input:
        shards = common_indel_sens_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = sens_prefix + ".common_sites.indels.bed.gz",
        mem_gb = 3.75,
        n_cpu = 2,
        docker = bcftools_docker
      }
  }

  # Collapse SV common site comparisons
  Array[File] common_sv_ppv_beds = select_all(CompareSvPpv.common_sites_bed)
  if ( length(common_sv_ppv_beds) > 0 ) {
    call QcTasks.ConcatTextFiles as CollapseCommonSvPpvSites {
      input:
        shards = common_sv_ppv_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = ppv_prefix + ".common_sites.svs.bed.gz",
        mem_gb = 3.75,
        n_cpu = 2,
        docker = bcftools_docker
      }
  }
  Array[File] common_sv_sens_beds = select_all(CompareSvSens.common_sites_bed)
  if ( length(common_sv_sens_beds) > 0 ) {
    call QcTasks.ConcatTextFiles as CollapseCommonSvSensSites {
      input:
        shards = common_sv_sens_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = sens_prefix + ".common_sites.svs.bed.gz",
        mem_gb = 3.75,
        n_cpu = 2,
        docker = bcftools_docker
      }
  }

  # Sum compressed concordance results by size
  call QcTasks.SumCompressedDistribs as SumPpvBySize {
    input:
      distrib_tsvs = flatten([select_all(CompareSnvPpv.val_by_size_tsv),
                              select_all(CompareIndelPpv.val_by_size_tsv),
                              select_all(CompareSvPpv.val_by_size_tsv)]),
      n_key_columns = 3,
      out_prefix = ppv_prefix + ".concordance_by_size",
      g2c_analysis_docker = g2c_analysis_docker
  }
  call QcTasks.SumCompressedDistribs as SumSensBySize {
    input:
      distrib_tsvs = flatten([select_all(CompareSnvSens.val_by_size_tsv),
                              select_all(CompareIndelSens.val_by_size_tsv),
                              select_all(CompareSvSens.val_by_size_tsv)]),
      n_key_columns = 3,
      out_prefix = sens_prefix + ".concordance_by_size",
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Sum compressed concordance results by frequency
  call QcTasks.SumCompressedDistribs as SumPpvByFreq {
    input:
      distrib_tsvs = flatten([select_all(CompareSnvPpv.val_by_freq_tsv),
                              select_all(CompareIndelPpv.val_by_freq_tsv),
                              select_all(CompareSvPpv.val_by_freq_tsv)]),
      n_key_columns = 3,
      out_prefix = ppv_prefix + ".concordance_by_af",
      g2c_analysis_docker = g2c_analysis_docker
  }
  call QcTasks.SumCompressedDistribs as SumSensByFreq {
    input:
      distrib_tsvs = flatten([select_all(CompareSnvSens.val_by_freq_tsv),
                              select_all(CompareIndelSens.val_by_freq_tsv),
                              select_all(CompareSvSens.val_by_freq_tsv)]),
      n_key_columns = 3,
      out_prefix = sens_prefix + ".concordance_by_af",
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Make full VID maps, if optioned
  if ( make_full_id_maps ) {
    # Collapse all site comparisons (not split by variant class)
    Array[File] all_ppv_shards = select_all(flatten([select_all(CompareSnvPpv.all_sites_bed),
                                                     select_all(CompareIndelPpv.all_sites_bed),
                                                     select_all(CompareSvPpv.all_sites_bed)]))
    call QcTasks.ConcatTextFiles as CollapseAllPpvSites {
      input:
        shards = all_ppv_shards,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = ppv_prefix + ".all_sites.bed.gz",
        mem_gb = 3.75,
        n_cpu = 2,
        docker = bcftools_docker
    }
    Array[File] all_sens_shards = select_all(flatten([select_all(CompareSnvSens.all_sites_bed),
                                                     select_all(CompareIndelSens.all_sites_bed),
                                                     select_all(CompareSvSens.all_sites_bed)]))
    call QcTasks.ConcatTextFiles as CollapseAllSensSites {
      input:
        shards = all_sens_shards,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = sens_prefix + ".all_sites.bed.gz",
        mem_gb = 3.75,
        n_cpu = 2,
        docker = bcftools_docker
    }

    # Extract consensus list mapping variant IDs between the two datasets
    call MakeVidMap as MakePpvVidMap {
      input:
        benchmark_bed = CollapseAllPpvSites.merged_file,
        out_prefix = ppv_prefix
    }
    call MakeVidMap as MakeSensVidMap {
      input:
        benchmark_bed = CollapseAllSensSites.merged_file,
        out_prefix = sens_prefix
    }
  }

  output {
    File? common_snv_ppv_bed = CollapseCommonSnvPpvSites.merged_file
    File? common_snv_sens_bed = CollapseCommonSnvSensSites.merged_file
    File? common_indel_ppv_bed = CollapseCommonIndelPpvSites.merged_file
    File? common_indel_sens_bed = CollapseCommonIndelSensSites.merged_file
    File? common_sv_ppv_bed = CollapseCommonSvPpvSites.merged_file
    File? common_sv_sens_bed = CollapseCommonSvSensSites.merged_file
    File ppv_by_size = SumPpvBySize.merged_distrib
    File sensitivity_by_size = SumSensBySize.merged_distrib
    File ppv_by_freq = SumPpvByFreq.merged_distrib
    File sensitivity_by_freq = SumSensByFreq.merged_distrib
    File? ppv_variant_id_map = MakePpvVidMap.vid_map
    File? sensitivity_variant_id_map = MakeSensVidMap.vid_map
  }
}


task CompareSites {
  input {
    File query_bed
    File ref_bed
    
    File genome_file
    String mode
    Float common_af
    Int overlap_pad = 1
    Boolean degenerate = false
    
    String prefix

    String? gcp_machine_type
    Int? disk_gb
    Float mem_gb = 4
    Int? n_cpu

    String g2c_analysis_docker
  }

  String machine_type = if defined(gcp_machine_type) then gcp_machine_type else ""
  Int n_cpu_default = floor(ceil(mem_gb) / 2)
  Int default_disk_gb = ceil(2 * size([query_bed, ref_bed], "GB")) + 20
  Int use_disk_gb = select_first([disk_gb, default_disk_gb])
  String degen_cmd = if degenerate then "--one-to-many" else ""

  command <<<
    set -eu -o pipefail

    # Start heartbeat to avoid silent VM death
    while true; do
      echo "[CompareSites] still running at $(date)"
      sleep 60
    done &
    HEARTBEAT_PID=$!
    trap "kill $HEARTBEAT_PID 2>/dev/null || true" EXIT

    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/compare_sites.py \
      -a ~{query_bed} \
      -b ~{ref_bed} \
      -g ~{genome_file} \
      -o ~{prefix} \
      --common-af ~{common_af} \
      --mode ~{mode} \
      --overlap-pad ~{overlap_pad} \
      ~{degen_cmd} \
      --no-reverse \
      --gzip
  >>>

  output {
    File all_sites_bed = "~{prefix}.loj.sites.bed.gz"
    File all_sites_bed_idx = "~{prefix}.loj.sites.bed.gz.tbi"
    File common_sites_bed = "~{prefix}.loj.sites.common.bed.gz"
    File common_sites_bed_idx = "~{prefix}.loj.sites.common.bed.gz.tbi"
    File val_by_size_tsv = "~{prefix}.loj.concordance_by_size.tsv.gz"
    File val_by_freq_tsv = "~{prefix}.loj.concordance_by_af.tsv.gz"
  }

  runtime {
    docker: g2c_analysis_docker
    predefinedMachineType: machine_type
    memory: "~{mem_gb} GB"
    cpu: select_first([n_cpu, n_cpu_default])
    disks: "local-disk ~{use_disk_gb} HDD"
    preemptible: 1
    max_retries: 1
  }
}


task MakeVidMap {
  input {
    File benchmark_bed
    String out_prefix
    Boolean invert = false

    Float mem_gb = 3.5
    String linux_docker = "ubuntu:plucky-20251001"
  }

  Int disk_gb = ceil(2 * size(benchmark_bed, "GB")) + 10
  String invert_cmd = if invert then "| awk -v OFS=\"\t\" '{ print $2, $1, $3 }' " else ""
  String out_fname = "~{out_prefix}.variant_id_map.tsv.gz"

  command <<<
    set -eu -o pipefail

    zcat ~{benchmark_bed} \
    | fgrep -v "#" \
    | awk -v FS="\t" \
      '{ if ($10=="NA") $11=100; printf "%s\t%s\t%f\n", $4, $9, $11 }' \
    | sort -nk3,3 \
    ~{invert_cmd} \
    | awk -v FS="\t" -v OFS="\t" '!seen[$1]++' \
    | gzip -c \
    > ~{out_fname} || true
  >>>

  output {
    File vid_map = out_fname
  }

  runtime {
    docker: linux_docker
    memory: "~{mem_gb} GB"
    cpu: 2
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}

