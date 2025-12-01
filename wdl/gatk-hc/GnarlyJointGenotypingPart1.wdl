# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# WDL to run first half of GATK "biggest practices" / "Gnarly" joint genotyping of SNVs/indels from gVCF input

# Unlike the parent WARP workflow, this task is designed to be executed on smaller subsets of sharded intervals
# E.g., per chromosome (or even smaller regions)
# These individual shards can then be merged prior to VQSR, fingerprinting, etc
# See GnarlyJointGenotypingPart2.wdl for those subsequent steps

# Based on GATK WARP pipeline:
# https://github.com/broadinstitute/warp/blob/develop/pipelines/broad/dna_seq/germline/joint_genotyping/JointGenotyping.wdl

# See also:
# https://gatk.broadinstitute.org/hc/en-us/articles/16957867036315-Introducing-GATK-Biggest-Practices-for-Joint-Calling-Supersized-Cohorts


version 1.0


import "JointGenotypingTasks.wdl" as Tasks
import "Utilities.wdl" as G2CUtils


workflow GnarlyJointGenotypingPart1 {
  input {
    File unpadded_intervals_file

    String callset_name
    File sample_name_map
    Int import_gvcfs_batch_size = 50

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File dbsnp_vcf

    Boolean make_hard_filtered_sites = true

    Int small_disk = 100
    Int import_gvcfs_disk_gb = 200
    Int make_sites_disk_gb = 200
    Int large_disk = 1000

    Int? top_level_scatter_count
    Boolean intervals_already_split = false
    Array[File]? custom_unpadded_intervals # If provided, this array will override any interval sharding operations
    Float unbounded_scatter_count_scale_factor = 0.15
    Int gnarly_scatter_count = 10
  }

  Array[Array[String]] sample_name_map_lines = read_tsv(sample_name_map)
  Int num_gvcfs = length(sample_name_map_lines)

  Int unbounded_scatter_count = select_first([top_level_scatter_count, round(unbounded_scatter_count_scale_factor * num_gvcfs)])
  Int scatter_count = if unbounded_scatter_count > 2 then unbounded_scatter_count else 2


  call Tasks.CheckSamplesUnique {
    input:
      sample_name_map = sample_name_map
  }

  # We deviate from the WARP protocol by enabling the user to provide custom pre-split intervals
  # This requires a custom task to strictly split the interval list into single-interval shards
  if ( ( intervals_already_split ) && ( !defined(custom_unpadded_intervals) ) ) {
    call G2CUtils.SplitIntervalList as G2CSplitIntervals {
      input:
        interval_list = unpadded_intervals_file
    }
  }
  if ( ( !intervals_already_split ) && ( !defined(custom_unpadded_intervals) ) ) {
    call Tasks.SplitIntervalList as GatkSplitIntervals {
      input:
        interval_list = unpadded_intervals_file,
        scatter_count = scatter_count,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        disk_size_gb = small_disk,
        sample_names_unique_done = CheckSamplesUnique.samples_unique
    }
  }
  
  Array[File] unpadded_intervals = select_first([custom_unpadded_intervals, 
                                                 G2CSplitIntervals.output_intervals, 
                                                 GatkSplitIntervals.output_intervals])

  scatter (idx in range(length(unpadded_intervals))) {
    call ImportGVCFsFT {
      input:
        sample_name_map = sample_name_map,
        interval = unpadded_intervals[idx],
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        workspace_dir_name = "genomicsdb",
        seed = idx,
        disk_size_gb = import_gvcfs_disk_gb,
        batch_size = import_gvcfs_batch_size
    }

    if ( gnarly_scatter_count > 1 ) {
      call Tasks.SplitIntervalList as GnarlyIntervalScatterDude {
        input:
          interval_list = unpadded_intervals[idx],
          scatter_count = gnarly_scatter_count,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          disk_size_gb = small_disk,
          sample_names_unique_done = CheckSamplesUnique.samples_unique
      }
    }

    Array[File] gnarly_intervals = select_first([GnarlyIntervalScatterDude.output_intervals, [unpadded_intervals[idx]]])

    scatter (gnarly_idx in range(length(gnarly_intervals))) {
      call GnarlyGenotyperFT {
        input:
          workspace_tar = ImportGVCFsFT.output_genomicsdb,
          interval = gnarly_intervals[gnarly_idx],
          output_vcf_filename = callset_name + "." + idx + "." + gnarly_idx + ".vcf.gz",
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          dbsnp_vcf = dbsnp_vcf
      }
    }

    Array[File] gnarly_vcfs = GnarlyGenotyperFT.output_vcf
    Array[File] gnarly_vcf_idxs = GnarlyGenotyperFT.output_vcf_index

    if ( length(gnarly_vcfs) > 1 ) {
      call Tasks.GatherVcfs as TotallyRadicalGatherVcfs {
        input:
          input_vcfs = gnarly_vcfs,
          output_vcf_name = callset_name + "." + idx + ".gnarly.vcf.gz",
          disk_size_gb = large_disk
      }
    }

    File genotyped_vcf = select_first([TotallyRadicalGatherVcfs.output_vcf, gnarly_vcfs[0]])
    File genotyped_vcf_index = select_first([TotallyRadicalGatherVcfs.output_vcf_index, gnarly_vcf_idxs[0]])

    if ( make_hard_filtered_sites ){
      call Tasks.HardFilterAndMakeSitesOnlyVcf {
        input:
          vcf = genotyped_vcf,
          vcf_index = genotyped_vcf_index,
          excess_het_threshold = 54.69,
          variant_filtered_vcf_filename = callset_name + "." + idx + ".variant_filtered.vcf.gz",
          sites_only_vcf_filename = callset_name + "." + idx + ".sites_only.variant_filtered.vcf.gz",
          disk_size_gb = make_sites_disk_gb
      }
    }

    File jg_vcf_shard = select_first([HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf, genotyped_vcf])
    File jg_vcf_shard_idx = select_first([HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index, genotyped_vcf_index])
  }

  output {
    Array[File] joint_genotyped_vcf = jg_vcf_shard
    Array[File] joint_genotyped_vcf_index = jg_vcf_shard_idx
    Array[File?] sites_only_vcfs = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf
    Array[File?] sites_only_vcfs_index = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf_index
  }
}


# Local version of ImportGVCFs to expose number of preemptible tries
# Also adds one max retry to account for Google API errors
# Also shuffles sample map order to improve performance when scalling to massive parallelization,
# as we have found that GCP I/O becomes glacially slow with too many parallel reads to the same URI
# Otherwise, this task is functionally identical to ImportGVCFs as imported in JointGenotypingTasks.wdl
# This task is suffixed with "FT" for "fault tolerant" to disambiguate
task ImportGVCFsFT {

  input {
    File sample_name_map
    File interval
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    String workspace_dir_name
    String seed

    Int disk_size_gb
    Int machine_mem_mb = 32000
    Int jvm_start_mb = 6000
    Int jvm_max_mb = 25000
    Int batch_size
    Int n_preemptible_tries = 1
    Int n_retries = 1

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
  }

  Int xms_mb_default = if machine_mem_mb / 3 < jvm_start_mb then jvm_start_mb else ceil(machine_mem_mb / 3)
  Int xms_mb = select_first([jvm_start_mb, xms_mb_default])
  Int xmx_mb_default = machine_mem_mb - 5000
  Int xmx_mb = select_first([jvm_max_mb, xmx_mb_default])
  String java_opts = "-Xms" + xms_mb + "m -Xmx" + xmx_mb + "m"

  command <<<
    set -euo pipefail

    rm -rf ~{workspace_dir_name}

    # Shuffle sample map to avoid I/O bottlenecking caused by
    # thousands of read operations on a single URI
    shuf --random-source=<( yes "~{seed}" ) ~{sample_name_map} > shuffled.sample.map.tsv

    # GATK dev comments below:
    # We've seen some GenomicsDB performance regressions related to intervals, so we're going to 
    # pretend we only have a single interval using the --merge-input-intervals arg
    # There's no data in between since we didn't run HaplotypeCaller over those loci 
    # so we're not wasting any compute

    # The memory setting here is very important and must be several GiB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    gatk --java-options "~{java_opts}" \
      GenomicsDBImport \
      --genomicsdb-workspace-path ~{workspace_dir_name} \
      --batch-size ~{batch_size} \
      -L ~{interval} \
      --sample-name-map shuffled.sample.map.tsv \
      --reader-threads 5 \
      --merge-input-intervals \
      --consolidate

    tar -cf ~{workspace_dir_name}.tar ~{workspace_dir_name}
  >>>

  runtime {
    memory: "~{machine_mem_mb} MiB"
    cpu: 4
    bootDiskSizeGb: 20
    disks: "local-disk " + disk_size_gb + " HDD"
    docker: gatk_docker
    preemptible: n_preemptible_tries
    max_retries: n_retries
  }

  output {
    File output_genomicsdb = "~{workspace_dir_name}.tar"
  }
}


# Local version of GnarlyGenotyper to expose number of preemptible tries
# Also adds one max retry to increase robustness to transiet errors
# Otherwise, this task is functionally identical to the one imported in JointGenotypingTasks.wdl
# This task is suffixed with "FT" for "fault tolerant" to disambiguate
task GnarlyGenotyperFT {

  input {
    File workspace_tar
    File interval
    String output_vcf_filename
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    String dbsnp_vcf
    Boolean make_annotation_db = false

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
    Int machine_mem_mb = 26000
    Int jvm_start_mb = 6000
    Int disk_size_gb = ceil(size(workspace_tar, "GiB") + size(ref_fasta, "GiB") + size(dbsnp_vcf, "GiB"))
    Int n_preemptible_tries = 3
  }

  Int xms_mb = if machine_mem_mb / 3 < jvm_start_mb then jvm_start_mb else ceil(machine_mem_mb / 3)
  Int xmx_mb = machine_mem_mb - 1000
  String java_opts = "-Xms" + xms_mb + "m -Xmx" + xmx_mb + "m"

  parameter_meta {
    interval: {
      localization_optional: true
    }
  }

  command <<<
    set -e

    tar -xf ~{workspace_tar}
    WORKSPACE=$( basename ~{workspace_tar} .tar)

    gatk --java-options "~{java_opts}" \
      GnarlyGenotyper \
      -R ~{ref_fasta} \
      -O ~{output_vcf_filename} \
      ~{true="--output-database-name annotationDB.vcf.gz" false="" make_annotation_db} \
      -D ~{dbsnp_vcf} \
      --only-output-calls-starting-in-intervals \
      -V gendb://$WORKSPACE \
      -L ~{interval} \
      -stand-call-conf 10 \
      --max-alternate-alleles 5 \
      --merge-input-intervals  
  >>>

  runtime {
    memory: "~{machine_mem_mb} MiB"
    cpu: 2
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size_gb + " HDD"
    preemptible: n_preemptible_tries
    max_retries: 1
    docker: gatk_docker
  }

  output {
    File output_vcf = "~{output_vcf_filename}"
    File output_vcf_index = "~{output_vcf_filename}.tbi"
    File? output_database = "annotationDB.vcf.gz"
    File? output_database_index = "annotationDB.vcf.gz.tbi"
  }
}

