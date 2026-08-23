# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Perform two-sided sample-level benchmarking of a "source" callset versus a reference "target" callset
# Note: this is a subroutine called within BenchmarkSamples.wdl, and is designed to be
# invoked by this parent workflow

# Expects that both callsets are provided according to G2C QC expectations


version 1.0


import "BenchmarkSites.wdl" as SiteBench
import "QcTasks.wdl" as QcTasks


workflow BenchmarkSamplesSingle {
  input {
    File source_all_sites_bed
    File? source_snv_bed
    File? source_indel_bed
    File? source_sv_bed
    File source_gt_tarball
    File source_sample_id_list
    String source_prefix

    File target_vcf
    File target_vcf_idx
    String target_prefix
    
    File id_map_tsv

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

  Int n_eval_intervals = length(eval_interval_beds)

  # Subset target VCF to overlapping samples
  call QcTasks.MakeHeaderFiller {}
  call SubsetTargetVcf {
    input:
      vcf = target_vcf,
      vcf_idx = target_vcf_idx,
      supp_header = MakeHeaderFiller.supp_vcf_header,
      id_map_tsv = id_map_tsv,
      source_sample_id_list = source_sample_id_list,
      bcftools_docker = bcftools_docker
  }

  # Only proceed if any overlapping samples were found
  if ( SubsetTargetVcf.n_overlapping_samples > 0 ) {

    # Collect site metrics for target VCF
    call QcTasks.CollectSiteMetrics as TargetVcfToBed {
      input:
        vcf = SubsetTargetVcf.filtered_sites_vcf,
        vcf_idx = SubsetTargetVcf.filtered_sites_vcf_idx,
        n_samples = SubsetTargetVcf.n_samples_in_original_vcf,
        common_af_cutoff = common_af_cutoff,
        g2c_analysis_docker = g2c_analysis_docker
    }

    # Run site-level benchmarking to determine mapping of 
    # matching variant IDs between source & target
    call SiteBench.BenchmarkSites {
      input:
        source_snv_bed = source_snv_bed,
        source_indel_bed = source_indel_bed,
        source_sv_bed = source_sv_bed,
        source_prefix = source_prefix,
        target_snv_bed = TargetVcfToBed.snv_sites,
        target_indel_bed = TargetVcfToBed.indel_sites,
        target_sv_bed = TargetVcfToBed.sv_sites,
        target_prefix = target_prefix,
        eval_interval_beds = eval_interval_beds,
        eval_interval_bed_names = eval_interval_bed_names,
        genome_file = genome_file,
        total_shards = total_shards,
        common_af_cutoff = common_af_cutoff,
        make_full_id_maps = true,
        snv_mem_scalar = snv_mem_scalar,
        indel_mem_scalar = indel_mem_scalar,
        bcftools_docker = bcftools_docker,
        g2c_analysis_docker = g2c_analysis_docker
    }

    # Collect sample genotypes for subsetted target VCF
    call QcTasks.CollectSampleGenotypeMetrics as GatherTargetGts {
      input:
        vcf = SubsetTargetVcf.filtered_dense_vcf,
        vcf_idx = SubsetTargetVcf.filtered_dense_vcf_idx,
        site_metrics = TargetVcfToBed.all_sites,
        g2c_analysis_docker = g2c_analysis_docker
    }
    call QcTasks.ConcatGenotypeTsvs as MakeTargetGtTarball {
      input:
        tsvs = [GatherTargetGts.genotypes_tsv],
        output_prefix = target_prefix,
        g2c_analysis_docker = g2c_analysis_docker
    }

    # Shard the cleaned ID map to optimally parallelize GT benchmarking tasks
    Int n_tasks_per_shard = 2 * n_eval_intervals
    Int n_naive_sample_splits = ceil(n_tasks_per_shard * SubsetTargetVcf.n_overlapping_samples / min_samples_per_shard)
    Int n_floored_sample_splits = if n_naive_sample_splits <= 1 then 1 else n_naive_sample_splits
    Int n_sample_shards = if n_naive_sample_splits < total_shards then n_floored_sample_splits else total_shards
    if ( n_sample_shards > 1 ) {
      call QcTasks.ShardTextFile as ShardSamples {
        input:
          input_file = SubsetTargetVcf.filtered_id_map,
          n_splits = n_sample_shards,
          out_prefix = "~{source_prefix}.~{target_prefix}.gt_bench.id_map_shard",
          shuffle = true,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
    Array[File] sample_shards = select_first([ShardSamples.shards, [SubsetTargetVcf.filtered_id_map]])

    # Keep benchmarking results separate for each evaluation interval set
    scatter ( i in range(n_eval_intervals) ) {

      String ei_prefix = if defined(eval_interval_bed_names)
                         then flatten(select_all([eval_interval_bed_names]))[i]
                         else "eval_interval_~{i}"
      String ppv_prefix = "~{ei_prefix}.~{source_prefix}.~{target_prefix}"
      String sens_prefix = "~{ei_prefix}.~{target_prefix}.~{source_prefix}"
      
      # Scatter over samples and perform benchmarking on each
      scatter ( k in range(length(sample_shards)) ) {

        # Benchmark PPV (source -> target)
        call QcTasks.BenchmarkGenotypes as BenchGtPpv {
          input:
            variant_id_map = select_all(BenchmarkSites.ppv_variant_id_maps)[i],
            sample_id_map = sample_shards[k],
            invert_sample_map = false,
            output_prefix = "~{ppv_prefix}.~{k}",
            source_gt_tarball = source_gt_tarball,
            target_gt_tarball = MakeTargetGtTarball.genotypes_tarball,
            source_site_metrics = source_all_sites_bed,
            common_af_cutoff = common_af_cutoff,
            g2c_analysis_docker = g2c_analysis_docker
        }

        # Benchmark sensitivity (target -> source)
        call QcTasks.BenchmarkGenotypes as BenchGtSens {
          input:
            variant_id_map = select_all(BenchmarkSites.sensitivity_variant_id_maps)[i],
            sample_id_map = sample_shards[k],
            invert_sample_map = true,
            output_prefix = "~{sens_prefix}.~{k}",
            source_gt_tarball = MakeTargetGtTarball.genotypes_tarball,
            target_gt_tarball = source_gt_tarball,
            source_site_metrics = TargetVcfToBed.all_sites,
            common_af_cutoff = common_af_cutoff,
            g2c_analysis_docker = g2c_analysis_docker
        }

      }

      # Collapse PPV benchmarking across all shards
      call QcTasks.SumCompressedDistribs as SumPpvDistrib {
        input:
          distrib_tsvs = BenchGtPpv.gt_bench_distrib,
          out_prefix = ppv_prefix + ".gt_comparison.distrib",
          n_key_columns = 5,
          g2c_analysis_docker = g2c_analysis_docker
      }

      # Collapse sensitivity benchmarking across all shards
      call QcTasks.SumCompressedDistribs as SumSensDistrib {
        input:
          distrib_tsvs = BenchGtSens.gt_bench_distrib,
          out_prefix = sens_prefix + ".gt_comparison.distrib",
          n_key_columns = 5,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
  }

  # If no overlapping samples are found, we still need to produce one empty 
  # file for each evalulation interval for each of PPV and sensitivity 
  # for downstream compatability
  if ( SubsetTargetVcf.n_overlapping_samples == 0 ) {

    # Keep benchmarking results separate for each evaluation interval set
    scatter ( i in range(n_eval_intervals) ) {

      String dummy_ei_prefix = if defined(eval_interval_bed_names)
                               then flatten(select_all([eval_interval_bed_names]))[i]
                               else "eval_interval_~{i}"
      String dummy_ppv_prefix = "~{dummy_ei_prefix}.~{source_prefix}.~{target_prefix}"
      String dummy_sens_prefix = "~{dummy_ei_prefix}.~{target_prefix}.~{source_prefix}"

      call MakeEmptyBenchDistrib as MakeDummyPpv {
        input:
          out_prefix = dummy_ppv_prefix + ".gt_comparison.distrib",
          bcftools_docker = bcftools_docker
      }

      call MakeEmptyBenchDistrib as MakeDummySens {
        input:
          out_prefix = dummy_sens_prefix + ".gt_comparison.distrib",
          bcftools_docker = bcftools_docker
      }
    }
  }


  output {
    # One file per evaluation interval set in each output array
    Array[File] ppv_distribs = select_first([SumPpvDistrib.merged_distrib, MakeDummyPpv.dummy_distrib])
    Array[File] sensitivity_distribs = select_first([SumSensDistrib.merged_distrib, MakeDummySens.dummy_distrib])
  }
}


# Subset a target VCF based on a dynamically-determined list of overlapping samples
task SubsetTargetVcf {
  input {
    File vcf
    File vcf_idx
    File supp_header
    File id_map_tsv
    File source_sample_id_list
    String bcftools_docker
    String gcp_machine_type = "n2d-standard-2"
  }

  Int disk_gb = ceil(2.5 * size(vcf, "GB")) + 10

  String filtered_map_fname = basename(id_map_tsv, ".tsv") + ".filtered.tsv"
  String filtered_dense_vcf_fname = basename(vcf, ".vcf.gz") + ".filtered.vcf.gz"
  String filtered_sites_vcf_fname = basename(vcf, ".vcf.gz") + ".filtered.sites.vcf.gz"

  command <<<
    set -eu -o pipefail

    # Get list of samples present in input VCF
    bcftools query -l ~{vcf} | sort > target.samples.list

    # Ensure all files are sorted for bash join
    sort ~{source_sample_id_list} > source.samples.list
    sort -k1,1 ~{id_map_tsv} > id_map.tsv

    # Use bash join to subset ID map
    join -j 1 -t $'\t' source.samples.list id_map.tsv \
    | sort -k2,2 \
    | join -1 2 -2 1 -t $'\t' - target.samples.list \
    | awk -v OFS="\t" '{ print $2, $1 }' \
    > ~{filtered_map_fname}

    # Subset input VCF to overlapping samples
    cut -f2 ~{filtered_map_fname} | sort -V | uniq > target.keep.list
    if [ $( cat target.keep.list | wc -l ) -lt 1 ]; then

      # If no overlapping samples exist, make empty VCFs for downstream compatability
      bcftools view \
        --header-only -G \
        -Oz -o ~{filtered_dense_vcf_fname} \
        ~{vcf}
      cp ~{filtered_dense_vcf_fname} ~{filtered_sites_vcf_fname}
      
    else

      # Make a dense VCF of overlapping samples
      bcftools annotate -h ~{supp_header} ~{vcf} \
      | bcftools view --no-update --samples-file target.keep.list --force-samples \
      | bcftools view --no-update --include 'GT="alt" | FILTER="MULTIALLELIC"' \
        -Oz -o ~{filtered_dense_vcf_fname}

      # Make a sites VCF for benchmarking comparison
      bcftools view --no-update -G \
        -Oz -o ~{filtered_sites_vcf_fname} \
        ~{filtered_dense_vcf_fname}
    
    fi

    # Index both versions of subsetted VCF
    tabix -p vcf -f ~{filtered_dense_vcf_fname}
    tabix -p vcf -f ~{filtered_sites_vcf_fname}
  >>>

  output {
    File filtered_dense_vcf = filtered_dense_vcf_fname
    File filtered_dense_vcf_idx = "~{filtered_dense_vcf_fname}.tbi"

    File filtered_sites_vcf = filtered_sites_vcf_fname
    File filtered_sites_vcf_idx = "~{filtered_sites_vcf_fname}.tbi"

    File filtered_id_map = filtered_map_fname

    Int n_samples_in_original_vcf = length(read_lines("target.samples.list"))
    Int n_overlapping_samples = length(read_lines("target.keep.list"))
  }

  runtime {
    docker: bcftools_docker
    predefinedMachineType: gcp_machine_type
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}


task MakeEmptyBenchDistrib {
  input {
    String out_prefix
    String bcftools_docker
  }

  command <<<
    set -eu -o pipefail

    echo -e "#sample\tclass\tsubclass\tfreq_bin\tzygosity\tno_match\tcarrier_match\tgt_match" \
    | gzip -c \
    > ~{out_prefix}.dummy.tsv.gz
  >>>

  output {
    File dummy_distrib = "~{out_prefix}.dummy.tsv.gz"
  }

  runtime {
    docker: bcftools_docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk 20 HDD"
    preemptible: 3
  }
}

