# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Unify GATK-HC and GATK-SV callsets while collapsing large indels & small SVs


version 1.0


import "NestedReshardVcfs.wdl" as NRV
import "Utilities.wdl" as Utils


workflow UnifyGatkCallsets {
  input {
    # Two ways to provide GATK-HC VCF information: as arrays for VCF & indexes, or
    # as a two-column .tsv with URIs for VCF and index. If both are provided,
    # the array-style inputs will be used.
    Array[File]? gatkhc_vcf_array
    Array[File]? gatkhc_vcf_idx_array
    File? gatkhc_vcf_info_tsv

    Array[File] gatksv_vcfs
    Array[File] gatksv_vcf_idxs

    File genome_file                               # BEDTools-style genome file
    Int min_interval_size = 1000000                # Minimum size of intervals used for SV/indel clustering
    Int intervals_per_shard_indel_clustering = 1   # Parallelization control for pre-clustering resharding task

    Int min_sv_size = 50
    Float size_scalar = 3                          # Maximum fold-difference between sizes of indels and SVs to tolerate
    String sv_mask_field = "SL"                    # Quality field to mask low-quality SV genotypes before determining matches to indels

    # BED4 files of final output intervals for resharding
    # Fourth column must correspond to desired output VCF name / interval name
    File snv_partition_intervals
    File indel_partition_intervals
    File sv_partition_intervals
    String large_sv_interval_name = "large_svs"    # Optional custom name for the extra "large SV" shard
    Int vcfs_per_shard = 20                        # Global parallelization control for *all* ReshardVcf tasks
    Int intervals_per_shard_final_partition = 3    # Parallelization control for final partitioning task
    Float final_partition_disk_scalar = 1.5        # Disk sizing parameter for final partitioning task

    String g2c_analysis_docker
    String linux_docker = "ubuntu:plucky-20251001"
  }

  Int min_large_indel_size = floor(min_sv_size / size_scalar)

  # Determine method of GATK-HC VCF input
  if ( !defined(gatkhc_vcf_array) || !defined(gatkhc_vcf_idx_array) ) {
    call Utils.ReadVcfInfo as ExtractGatkhcInfo {
      input:
        vcf_info = select_first([gatkhc_vcf_info_tsv, '']),
        linux_docker = linux_docker
    }
  }
  Array[File] gatkhc_vcfs = select_first([ExtractGatkhcInfo.vcf_uris, gatkhc_vcf_array])
  Array[File] gatkhc_vcf_idxs = select_first([ExtractGatkhcInfo.vcf_tbi_uris, gatkhc_vcf_idx_array])

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
  call Utils.Max as CalcMaxIndelSize {input: values = select_all(SplitGatkHcBySize.largest_indel_size)}

  # Pre-filter SVs to only those eligible for possible indel clustering
  Array[Pair[File, File]] sv_vcf_infos = zip(gatksv_vcfs, gatksv_vcf_idxs)
  scatter ( sv_vcf_info in sv_vcf_infos ) {
    call PrepareSvs {
      input:
        vcf = sv_vcf_info.left,
        vcf_idx = sv_vcf_info.right,
        max_size = ceil(size_scalar * CalcMaxIndelSize.max) + ceil(size_scalar),
        g2c_analysis_docker = g2c_analysis_docker
    }
  }

  # Make indel/sv clustering intervals
  call MakeClusteringIntervals {
    input:
      genome_file = genome_file,
      indel_beds = select_all(SplitGatkHcBySize.large_indel_sites_bed),
      sv_beds = select_all(PrepareSvs.eligible_sv_sites_bed),
      min_interval_size = min_interval_size,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Partition large indels across clustering intervals
  call NRV.NestedReshardVcfs as SplitIndelsForClustering {
    input:
      vcfs = select_all(SplitGatkHcBySize.large_indel_vcf),
      vcf_idxs = select_all(SplitGatkHcBySize.large_indel_vcf_idx),
      reshard_intervals_bed = MakeClusteringIntervals.intervals_bed,
      interval_suffix = "large_indels",
      rename_variants = true,
      vcfs_per_shard = vcfs_per_shard,
      intervals_per_shard = intervals_per_shard_indel_clustering,
      reshard_task_mem_gb = 12,
      reshard_task_n_cpu = 4,
      g2c_analysis_docker = g2c_analysis_docker,
      linux_docker = linux_docker
  }

  # Partition SVs across clustering intervals
  call Utils.ReshardVcfs as SplitSvsForClustering {
    input:
      vcfs = PrepareSvs.eligible_sv_vcf,
      vcf_idxs = PrepareSvs.eligible_sv_vcf_idx,
      intervals_bed = MakeClusteringIntervals.intervals_bed,
      interval_suffix = "svs",
      rename = true,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Sort partitioned indels and SVs to ensure no mismatching due to wildcard globs
  call AlignIndelAndSvFiles {
    input:
      indel_vcfs = SplitIndelsForClustering.resharded_vcfs,
      indel_vcf_idxs = SplitIndelsForClustering.resharded_vcf_idxs,
      sv_vcfs = select_all(SplitSvsForClustering.resharded_vcfs),
      sv_vcf_idxs = select_all(SplitSvsForClustering.resharded_vcf_idxs),
      interval_names = MakeClusteringIntervals.interval_names_list,
      linux_docker = linux_docker
  }

  # Cluster large indels and SVs
  scatter ( interval_info in AlignIndelAndSvFiles.interval_infos ) {

    call DefineClusters {
      input:
        indel_vcf = interval_info[1],
        indel_vcf_idx = interval_info[2],
        sv_vcf = interval_info[3],
        sv_vcf_idx = interval_info[4],
        sv_mask_field = sv_mask_field,
        genome_file = genome_file,
        size_scalar = size_scalar,
        g2c_analysis_docker = g2c_analysis_docker
    }

    call ResolveClusters {
      input:
        indel_vcf = interval_info[1],
        indel_vcf_idx = interval_info[2],
        sv_vcf = interval_info[3],
        sv_vcf_idx = interval_info[4],
        clusters = DefineClusters.final_clusters,
        combined_header = DefineClusters.combined_header,
        g2c_analysis_docker = g2c_analysis_docker
    }

  }
  call Utils.ConcatTextFiles as ConcatClusterLogs {
    input:
      shards = DefineClusters.final_clusters,
      concat_command = "zcat",
      compression_command = "gzip -c",
      input_has_header = true,
      output_filename = "indel_sv_cluster_assignments.tsv.gz",
      docker = linux_docker
  }

  # Reheader small indels to be compatible with reclustered large indels
  call Utils.GetVcfHeader as GetFinalIndelHeader {
    input:
      vcf = select_first(select_all(ResolveClusters.integrated_indel_vcf)),
      vcf_idx = select_first(select_all(ResolveClusters.integrated_indel_vcf_idx)),
      bcftools_docker = g2c_analysis_docker
  }
  Array[Pair[File, File]] si_infos = zip(select_all(SplitGatkHcBySize.small_indel_vcf),
                                         select_all(SplitGatkHcBySize.small_indel_vcf_idx))
  scatter ( si_info in si_infos ) {
    call Utils.ReheaderVcf as ReheaderSmallIndels {
      input:
        vcf = si_info.left,
        vcf_idx = si_info.right,
        new_header = GetFinalIndelHeader.header,
        bcftools_docker = g2c_analysis_docker
    }
  }

  # Reheader passthrough SVs to be compatible with reclustered small SVs
  call Utils.GetVcfHeader as GetFinalSvHeader {
    input:
      vcf = select_first(select_all(ResolveClusters.integrated_sv_vcf)),
      vcf_idx = select_first(select_all(ResolveClusters.integrated_sv_vcf_idx)),
      bcftools_docker = g2c_analysis_docker
  }
  Array[Pair[File, File]] pt_infos = zip(PrepareSvs.passthrough_sv_vcf,
                                         PrepareSvs.passthrough_sv_vcf_idx)
  scatter ( pt_info in pt_infos ) {
    call Utils.ReheaderVcf as ReheaderPassthroughSvs {
      input:
        vcf = pt_info.left,
        vcf_idx = pt_info.right,
        new_header = GetFinalSvHeader.header,
        bcftools_docker = g2c_analysis_docker
    }
  }

  # Postprocess SNVs
  call NRV.NestedReshardVcfs as PartitionSnvOutputs {
    input:
      vcfs = select_all(SplitGatkHcBySize.snv_vcf),
      vcf_idxs = select_all(SplitGatkHcBySize.snv_vcf_idx),
      reshard_intervals_bed = snv_partition_intervals,
      rename_variants = true,
      keep_empty_resharded_vcfs = false,
      vcfs_per_shard = vcfs_per_shard,
      intervals_per_shard = intervals_per_shard_final_partition,
      reshard_task_mem_gb = 12,
      reshard_task_n_cpu = 4,
      g2c_analysis_docker = g2c_analysis_docker,
      linux_docker = linux_docker
  }

  # Postprocess indels
  Array[File] all_indel_vcfs = select_all(flatten([ReheaderSmallIndels.reheadered_vcf,
                                                   select_all(ResolveClusters.integrated_indel_vcf)]))
  Array[File] all_indel_vcf_idxs = select_all(flatten([ReheaderSmallIndels.reheadered_vcf_idx,
                                                       select_all(ResolveClusters.integrated_indel_vcf_idx)]))
  call NRV.NestedReshardVcfs as PartitionIndelOutputs {
    input:
      vcfs = all_indel_vcfs,
      vcf_idxs = all_indel_vcf_idxs,
      reshard_intervals_bed = indel_partition_intervals,
      resharded_vcf_header = select_first(select_all(ResolveClusters.integrated_indel_vcf)),
      rename_variants = true,
      keep_empty_resharded_vcfs = false,
      vcfs_per_shard = vcfs_per_shard,
      intervals_per_shard = intervals_per_shard_final_partition,
      reshard_task_mem_gb = 12,
      reshard_task_n_cpu = 4,
      g2c_analysis_docker = g2c_analysis_docker,
      linux_docker = linux_docker
  }

  # Postprocess SVs
  Array[File] all_sv_vcfs = select_all(flatten([ReheaderPassthroughSvs.reheadered_vcf,
                                                ResolveClusters.integrated_sv_vcf]))
  Array[File] all_sv_vcf_idxs = select_all(flatten([ReheaderPassthroughSvs.reheadered_vcf_idx,
                                                    ResolveClusters.integrated_sv_vcf_idx]))
  call NRV.NestedReshardVcfs as PartitionSvOutputs {
    input:
      vcfs = all_sv_vcfs,
      vcf_idxs = all_sv_vcf_idxs,
      reshard_intervals_bed = sv_partition_intervals,
      resharded_vcf_header = select_first(select_all(ResolveClusters.integrated_sv_vcf)),
      rename_variants = true,
      keep_empty_resharded_vcfs = false,
      vcfs_per_shard = vcfs_per_shard,
      intervals_per_shard = intervals_per_shard_final_partition,
      g2c_analysis_docker = g2c_analysis_docker,
      linux_docker = linux_docker
  }

  # Make one separate SV shard for all large SVs to aid in downstream region-based queries
  call DefineLargeSVCutoff {
    input:
      intervals = sv_partition_intervals,
      docker = linux_docker
  }
  Array[File] partitioned_sv_vcfs = PartitionSvOutputs.resharded_vcfs
  Array[File] partitioned_sv_vcf_idxs = PartitionSvOutputs.resharded_vcf_idxs
  Array[Pair[File, File]] penultimate_sv_infos = zip(partitioned_sv_vcfs, partitioned_sv_vcf_idxs)
  scatter ( svpair in penultimate_sv_infos ) {
    call ExtractLargeSvs {
      input:
        vcf = svpair.left,
        vcf_idx = svpair.right,
        size_cutoff = DefineLargeSVCutoff.cutoff,
        bcftools_docker = g2c_analysis_docker
    }
  }
  Array[File] large_sv_vcfs = select_all(ExtractLargeSvs.large_vcf)
  Array[File] large_sv_vcf_idxs = select_all(ExtractLargeSvs.large_vcf_idx)
  if ( length(large_sv_vcfs) > 0 ) {
    call Utils.ConcatVcfs as ConcatLargeSvs {
      input:
        vcfs = large_sv_vcfs,
        vcf_idxs = large_sv_vcf_idxs,
        out_prefix = large_sv_interval_name,
        bcftools_concat_options = "-a",
        bcftools_docker = g2c_analysis_docker
    }
  }

  # Prepare reporting log of all variants considered & collapsed
  call Utils.Sum as SumSnvCounts {input: values = select_all(SplitGatkHcBySize.n_snvs)}
  call Utils.Sum as SumSmallIndelCounts {input: values = select_all(SplitGatkHcBySize.n_small_indels)}
  call Utils.Sum as SumLargeIndelCounts {input: values = select_all(SplitGatkHcBySize.n_large_indels)}
  call Utils.Sum as SumLargeIndelGE50Counts {input: values = select_all(SplitGatkHcBySize.n_large_indels_ge50bp)}
  call Utils.Sum as SumEligSvCounts {input: values = select_all(PrepareSvs.n_eligible_svs)}
  call Utils.Sum as SumPtSvCounts {input: values = select_all(PrepareSvs.n_passthrough_svs)}
  call Utils.Sum as SumClusterCounts {input: values = select_all(DefineClusters.n_clusters)}
  call Utils.Sum as SumMultiIndelClusterCounts {input: values = select_all(DefineClusters.n_multi_indel_clusters)}
  call Utils.Sum as SumClusteredIndelCounts {input: values = select_all(DefineClusters.n_clustered_indels)}
  call Utils.Sum as SumMultiSvClusterCounts {input: values = select_all(DefineClusters.n_multi_sv_clusters)}
  call Utils.Sum as SumClusteredSvCounts {input: values = select_all(DefineClusters.n_clustered_svs)}
  call WriteSummaryLog {
    input:
      snv_count = ceil(SumSnvCounts.sum),
      small_indel_count = ceil(SumSmallIndelCounts.sum),
      large_indel_count = ceil(SumLargeIndelCounts.sum),
      large_indel_ge50bp_count = ceil(SumLargeIndelGE50Counts.sum),
      eligible_sv_count = ceil(SumEligSvCounts.sum),
      passthrough_sv_count = ceil(SumPtSvCounts.sum),
      cluster_count = ceil(SumClusterCounts.sum),
      multi_indel_cluster_count = ceil(SumMultiIndelClusterCounts.sum),
      clustered_indel_count = ceil(SumClusteredIndelCounts.sum),
      multi_sv_cluster_count = ceil(SumMultiSvClusterCounts.sum),
      clustered_sv_count = ceil(SumClusteredSvCounts.sum),
      linux_docker = linux_docker
  }

  output {
    Array[File?] cleaned_snv_vcfs = PartitionSnvOutputs.resharded_vcfs
    Array[File?] cleaned_snv_vcf_idxs = PartitionSnvOutputs.resharded_vcf_idxs

    Array[File?] cleaned_indel_vcfs = PartitionIndelOutputs.resharded_vcfs
    Array[File?] cleaned_indel_vcf_idxs = PartitionIndelOutputs.resharded_vcf_idxs

    Array[File?] cleaned_sv_vcfs = select_all(flatten([ExtractLargeSvs.notlarge_vcf, [ConcatLargeSvs.merged_vcf]]))
    Array[File?] cleaned_sv_vcf_idxs = select_all(flatten([ExtractLargeSvs.notlarge_vcf_idx, [ConcatLargeSvs.merged_vcf_idx]]))

    File cluster_assignment_log = ConcatClusterLogs.merged_file
    File integration_summary_log = WriteSummaryLog.logfile
  }
}


# Make intervals for indel/SV clustering based on N-masked reference regions
task MakeClusteringIntervals {
  input {
    File genome_file
    Array[File] indel_beds
    Array[File] sv_beds
    Int min_interval_size
    Int min_gap_size = 500
    String g2c_analysis_docker
  }

  Array[File] all_beds = flatten([indel_beds, sv_beds])
  Int buffer_width = floor(min_gap_size / 2)

  command <<<
    set -eu -o pipefail

    # Merge indel and SV coordinates
    cat ~{write_lines(all_beds)} \
    | xargs -I {} zcat {} \
    | grep -ve '^#' \
    | awk -v buf="~{buffer_width}" -v FS="\t" -v OFS="\t" \
      '{pos=$2-buf; if (pos<0) pos=0; print $1, pos, $3+buf }' \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | bedtools merge -i - \
    | bgzip -c \
    > variant_coverage.bed.gz

    # Define eligible split points
    awk -v FS="\t" -v OFS="\t" '{ print $1, 0, $2 }' ~{genome_file} \
    | bedtools subtract -a - -b variant_coverage.bed.gz \
    | awk -v ms="~{min_gap_size}" -v FS="\t" -v OFS="\t" \
      '{ if ($3-$2>=ms) print }' \
    | bgzip -c \
    > eligible_splitpoints.bed.gz

    # Filter splitpoints to retain only those at least $min_interval_size apart
    zcat eligible_splitpoints.bed.gz \
    | awk -v md="~{min_interval_size}" -v FS="\t" -v OFS="\t" \
      'NR==1 { prev_chr=$1; prev_end=$3; print; next }
      { if ($1 != prev_chr || $2 >= prev_end + md) { prev_chr=$1; prev_end=$3; print } }' \
    | bgzip -c \
    > splits.filtered.bed.gz

    # Subtract qualifying splitpoints from genome file to define intervals
    awk -v FS="\t" -v OFS="\t" '{ print $1, 0, $2 }' ~{genome_file} \
    | bedtools subtract -a - -b splits.filtered.bed.gz \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | awk -v OFS="\t" '{ print $0, "clustering_interval_"NR }' \
    | bgzip -c \
    > clustering.intervals.bed.gz

    # Write list of clustering interval names to be used later
    zcat clustering.intervals.bed.gz \
    | cut -f4 > "clustering.intervals.ids.list"
  >>>

  output {
    File intervals_bed = "clustering.intervals.bed.gz"
    File interval_names_list = "clustering.intervals.ids.list"
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
    Boolean strict = true
    Float mem_gb = 1.75
    Int n_cpu = 1
    String g2c_analysis_docker
  }

  Int disk_gb = ceil(3 * size(vcf, "GB")) + 15
  String out_prefix = basename(vcf, ".vcf.gz")
  String snv_out_vcf = out_prefix + ".snv.vcf.gz"
  String si_out_vcf = out_prefix + ".small_indel.vcf.gz"
  String li_out_vcf = out_prefix + ".large_indel.vcf.gz"
  String li_out_sites_bed = out_prefix + ".large_indel.sites.bed.gz"

  command <<<
    set -eu -o pipefail

    # Ensure tabix index is localized to the same directory as the input VCF
    if [ "~{vcf}.tbi" != "~{vcf_idx}" ]; then
      mv "~{vcf_idx}" "~{vcf}.tbi"
    fi

    # Initialize trackers
    for fname in snv.variant_count.txt \
                 small_indel.variant_count.txt \
                 large_indel.variant_count.txt \
                 large_indel.ge50bp.variant_count.txt \
                 "~{out_prefix}.largest_indel_size.txt"; do
      echo "0" > $fname
    done

    # If input VCF is empty, we can just exit early
    if [ $( bcftools index -n ~{vcf} ) -eq 0 ]; then
      exit 0
    fi

    # Split input VCF
    /opt/pancan_germline_wgs/scripts/gatkhc_helpers/split_gatkhc_by_size.py \
      -i ~{vcf} \
      --large-indel-min-size ~{large_indel_size} \
      --largest-indel-log "~{out_prefix}.largest_indel_size.txt" \
      --make-large-indel-bed \
      -o ~{out_prefix}
    if [ -s "~{out_prefix}.large_indel.sites.bed" ]; then
      sort -Vk1,1 -k2,2n -k3,3n "~{out_prefix}.large_indel.sites.bed" \
      | cat <( echo -e "#chrom\tstart\tend" ) - \
      | bgzip -c \
      > "~{li_out_sites_bed}"
      zcat "~{li_out_sites_bed}" \
      | fgrep -v "#" \
      | awk -v FS="\t" '{ if ($3-$2>=50) print }' \
      | wc -l \
      > "large_indel.ge50bp.variant_count.txt" || true
    fi

    # Index and count variants for each output VCF, and delete empty VCFs
    for vc in snv small_indel large_indel; do
      vk=0
      ovcf="~{out_prefix}.$vc.vcf.gz"
      if [ -e $ovcf ]; then
        tabix -p vcf -f $ovcf
        vk=$( bcftools index -n $ovcf )
      fi
      if [ $vk -eq 0 ]; then
        rm $ovcf $ovcf.tbi
      fi
      echo $vk > $vc.variant_count.txt
    done

    # Validation check: since all outputs are optional, this task can seemingly
    # be interpreted by Cromwell as a success even if it errors out. Thus,
    # we need to explicitly exit non-zero if there are no VCFs in pwd.
    # This should be safe since we ruled out empty input VCFs earlier
    if [ "~{strict}" == "true" ]; then
      if [ $( find ./ -name "~{out_prefix}.*.vcf.gz" | wc -l ) -lt 1 ]; then
        echo "Likely error: no output VCFs were generated. This is usually unexpected and indicates a silent error."
        sleep 60s
        exit 1
      fi
    fi

    # Sometimes this task seems to finish unusually quickly, and Cromwell
    # seemingly doesn't always catch this in highly parallelized use-cases.
    # To be safe, we can sleep for 30 seconds to allow Cromwell to catch up
    sleep 30s
  >>>

  output {
    File? snv_vcf = snv_out_vcf
    File? snv_vcf_idx = snv_out_vcf + ".tbi"
    Int? n_snvs = read_int("snv.variant_count.txt")

    File? small_indel_vcf = si_out_vcf
    File? small_indel_vcf_idx = si_out_vcf + ".tbi"
    Int? n_small_indels = read_int("small_indel.variant_count.txt")

    File? large_indel_vcf = li_out_vcf
    File? large_indel_vcf_idx = li_out_vcf + ".tbi"
    File? large_indel_sites_bed = li_out_sites_bed
    Int? n_large_indels = read_int("large_indel.variant_count.txt")
    Int? n_large_indels_ge50bp = read_int("large_indel.ge50bp.variant_count.txt")

    Int largest_indel_size = read_int("~{out_prefix}.largest_indel_size.txt")
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


task PrepareSvs {
  input {
    File vcf
    File vcf_idx
    Int max_size
    String g2c_analysis_docker
  }

  String output_prefix = basename(vcf, ".vcf.gz")
  String elig_outfile = output_prefix + ".eligible_svs.vcf.gz"
  String elig_sites_bed = output_prefix + ".eligible_svs.sites.bed.gz"
  String pt_outfile = output_prefix + ".passthrough_svs.vcf.gz"
  Int disk_gb = (3 * ceil(size([vcf], "GB"))) + 10

  command <<<
    set -eu -o pipefail

    # Although slightly less efficient, this is clean enough to apply as a double-pass with bcftools

    # Inclusion
    bcftools view \
      -i 'INFO/SVTYPE = "DEL,DUP,INS" & INFO/SVLEN <= ~{max_size}' \
      -Oz -o "~{elig_outfile}" \
      ~{vcf}
    tabix -p vcf "~{elig_outfile}"
    bcftools index -n "~{elig_outfile}" > eligible.count.txt

    # Get site coordinates for eligible SVs to help with interval sharding downstream
    if [ $( bcftools index -n "~{elig_outfile}" ) -gt 0 ]; then
      bcftools query \
        -f '%CHROM\t%POS\t%END\t%REF\t%ALT\t%INFO/SVLEN\t.\t.\t.\t.\t.\t.\t.\t.\t.\n' \
         "~{elig_outfile}" \
      | /opt/pancan_germline_wgs/scripts/qc/vcf_qc/clean_site_metrics.py \
        -o eligible --gzip -N 100
      zcat eligible.sv.sites.bed.gz \
      | cut -f1-3 | bgzip -c \
      > "~{elig_sites_bed}"
    fi

    # Exclusion
    bcftools view \
      -e 'INFO/SVTYPE = "DEL,DUP,INS" & INFO/SVLEN <= ~{max_size}' \
      -Oz -o "~{pt_outfile}" \
      ~{vcf}
    tabix -p vcf "~{pt_outfile}"
    bcftools index -n "~{pt_outfile}" > passthrough.count.txt
  >>>

  output {
    File eligible_sv_vcf = "~{elig_outfile}"
    File eligible_sv_vcf_idx = "~{elig_outfile}.tbi"
    File? eligible_sv_sites_bed = "~{elig_sites_bed}"
    Int n_eligible_svs = read_int("eligible.count.txt")

    File passthrough_sv_vcf = "~{pt_outfile}"
    File passthrough_sv_vcf_idx = "~{pt_outfile}.tbi"
    Float passthrough_sv_vcf_size = size(passthrough_sv_vcf, "GB")
    Int n_passthrough_svs = read_int("passthrough.count.txt")
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 1
    maxRetries: 1
  }
}


# Ensure sharded indel and SV files are sorted identically for subsequent clustering
task AlignIndelAndSvFiles {
  input {
    Array[String] indel_vcfs
    Array[String] indel_vcf_idxs
    Array[String] sv_vcfs
    Array[String] sv_vcf_idxs
    File interval_names
    String linux_docker
  }

  command <<<
    set -eu -o pipefail

      cat ~{write_lines(indel_vcfs)} > indel_vcf.uris.list
      cat ~{write_lines(indel_vcf_idxs)} > indel_tbi.uris.list
      cat ~{write_lines(sv_vcfs)} > sv_vcf.uris.list
      cat ~{write_lines(sv_vcf_idxs)} > sv_tbi.uris.list

    while read iid; do
      for wrapper in 1; do
        echo $iid
        fgrep -w "$iid.large_indels.sorted.vcf.gz" indel_vcf.uris.list
        fgrep -w "$iid.large_indels.sorted.vcf.gz.tbi" indel_tbi.uris.list
        fgrep -w "$iid.svs.vcf.gz" sv_vcf.uris.list
        fgrep -w "$iid.svs.vcf.gz.tbi" sv_tbi.uris.list
      done | paste - - - - - >> ordered_interval_info.tsv
    done < ~{interval_names}
  >>>

  output {
    Array[Array[String]] interval_infos = read_tsv("ordered_interval_info.tsv")
  }

  runtime {
    docker: linux_docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk 20 HDD"
    preemptible: 1
    maxRetries: 3
  }
}


# Define clusters of indels and SVs to be integrated
task DefineClusters {
  input {
    File indel_vcf
    File indel_vcf_idx
    
    File sv_vcf
    File sv_vcf_idx
    String sv_mask_field = "SL"
    
    File genome_file
    Float size_scalar
    Float min_nonref_jaccard = 0.1

    String g2c_analysis_docker

    Float mem_gb = 3.5
    Int n_cpu = 2
  }

  Int disk_gb = ceil(2.5 * size([sv_vcf, indel_vcf], "GB")) + 20
  String output_prefix = basename(sv_vcf, ".svs.vcf.gz")
  String out_fname = "~{output_prefix}.final_clusters.tsv.gz"

  command <<<
    set -eu -o pipefail

    # Start heartbeat to avoid silent VM death
    (
      while true; do
        echo "[DefineClusters] still running at $(date)"
        sleep 60
      done
    ) &
    HEARTBEAT_PID=$!
    trap "kill $HEARTBEAT_PID 2>/dev/null || true" EXIT

    # Initialize trackers
    for fname in total_clusters.count.txt \
                multi_indel.clusters.count.txt \
                indel.variant_count.txt \
                multi_sv.clusters.count.txt \
                sv.variant_count.txt; do
      echo "0" > $fname
    done

    # Ensure VCF indexes localize to same directory as VCFs
    if [ "~{indel_vcf}.tbi" != "~{indel_vcf_idx}" ]; then
      cp "~{indel_vcf_idx}" "~{indel_vcf}.tbi"
    fi
    if [ "~{sv_vcf}.tbi" != "~{sv_vcf_idx}" ]; then
      cp "~{sv_vcf_idx}" "~{sv_vcf}.tbi"
    fi

    # Make header as helper for resolution step
    bcftools query -l ~{indel_vcf} > indel.samples.list
    bcftools query -l ~{sv_vcf} > sv.samples.list
    cat indel.samples.list sv.samples.list | sort -V | uniq > union.samples.list
    bcftools view \
      --samples-file union.samples.list \
      --force-samples \
      --header-only \
      -Oz -o indel.header.vcf.gz \
      ~{indel_vcf}
    bcftools view \
      --samples-file union.samples.list \
      --force-samples \
      --header-only \
      -Oz -o sv.header.vcf.gz \
      ~{sv_vcf}
    bcftools concat \
      indel.header.vcf.gz \
      sv.header.vcf.gz \
    | bcftools view \
      --samples-file union.samples.list \
      --force-samples \
      --header-only \
      -Oz -o combined.header.vcf.gz
    tabix -p vcf -f combined.header.vcf.gz

    # Step 0: if zero indels or SVs are present, there's nothing to do other than make dummy files
    n_indel=$( bcftools index -n ~{indel_vcf} )
    n_sv=$( bcftools index -n ~{sv_vcf} )
    if [ $n_indel -eq 0 ] || [ $n_sv -eq 0 ]; then
      echo -e "#SVs\tindels" | gzip -c > "~{out_fname}"
      exit 0
    fi

    # Step 1: Extract BED info for indels and SVs
    echo "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">" > header.supp.vcf
    for vc in indel sv; do
      case "$vc" in
        "indel")
          invcf="~{indel_vcf}"
          ;;
        "sv")
          invcf="~{sv_vcf}"
          ;;
      esac
      bcftools view $invcf \
      | bcftools annotate -h header.supp.vcf \
      | bcftools +fill-tags -- -t AC,AN,AF \
      | bcftools query \
        -f '%CHROM\t%POS\t%END\t%REF\t%ALT\t%INFO/SVLEN\t%INFO/AN\t%INFO/AC\t%INFO/AF\t.\t.\t.\t.\t.\t.\n' \
      | /opt/pancan_germline_wgs/scripts/qc/vcf_qc/clean_site_metrics.py \
        -o input.$vc \
        --gzip \
        -N 100
      find ./ -name "input.$vc.*.sites.bed.gz" > $vc.partitioned.sites.list
      n_parts=$( wc -l < $vc.partitioned.sites.list )
      if [ $n_parts -eq 0 ]; then
        echo "ERROR: no partitioned BEDs found for $vc" >&2
        exit 1
      elif [ $n_parts -gt 1 ]; then
        cat \
          <( zcat $( sed -n '1p' $vc.partitioned.sites.list ) | sed -n '1p' ) \
          <( cat $vc.partitioned.sites.list | xargs -I {} zcat {} \
             | grep -ve '^#' | sort -Vk1,1 -k2,2n -k3,3n -k4,4V ) \
        | bgzip -c \
        > input.$vc.sites.bed.gz
      else
        cp "$( sed -n '1p' $vc.partitioned.sites.list )" input.$vc.sites.bed.gz
      fi
      tabix -p bed -f input.$vc.sites.bed.gz
    done

    # Step 2a: find candidate overlapping indels & SVs
    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/compare_sites.py \
      -a input.sv.sites.bed.gz \
      -b input.indel.sites.bed.gz \
      -g ~{genome_file} \
      -o candidate_hits \
      --mode both \
      --max-overlap-size-diff 3.0 \
      --overlap-pad 1 \
      --one-to-many \
      --gzip
    zcat candidate_hits.loj.sites.bed.gz \
    | grep -ve '^#' \
    | awk -v FS="\t" -v OFS="\t" '{ if ($9!="NA") print "sv_"$4, "indel_"$9, $NF }' \
    > sv_to_indel.candidates.tsv
    zcat candidate_hits.roj.sites.bed.gz \
    | grep -ve '^#' \
    | awk -v FS="\t" -v OFS="\t" '{ if ($9!="NA") print "sv_"$9, "indel_"$4, $NF }' \
    > indel_to_sv.candidates.tsv
    cat indel_to_sv.candidates.tsv sv_to_indel.candidates.tsv \
    | sort -Vk1,1 -k3,3nr -k2,2V \
    | uniq \
    | gzip -c \
    > candidate_hits.pairs.tsv.gz
    zcat candidate_hits.pairs.tsv.gz \
    | awk -v OFS="\n" '{ print $1, $2 }' \
    | sort -V \
    | uniq \
    > candidate_hits.vids.list

    # Step 2b: if no candidate hits are found, there's nothing left to do other than make dummy files
    n_candidates=$( wc -l < candidate_hits.vids.list )
    if [ $n_candidates -eq 0 ]; then
      echo -e "#SVs\tindels" | gzip -c > "~{out_fname}"
      exit 0
    fi

    # Step 3: Compute non-ref GT Jaccard indexes for all candidate pairs
    # Note that the double call of SV GT masking is *not* by mistake; it is intentional
    bcftools view \
      --samples-file union.samples.list \
      --force-samples \
      "~{indel_vcf}" \
    | bcftools annotate --set-id "indel_%ID" \
    | bcftools view \
      --include 'ID=@candidate_hits.vids.list' \
      -Oz -o indel.for_jaccard.vcf.gz
    bcftools view \
      --samples-file union.samples.list \
      --force-samples \
      "~{sv_vcf}" \
    | bcftools annotate --set-id "sv_%ID" \
    | bcftools view \
      --include 'ID=@candidate_hits.vids.list' \
    | /opt/pancan_germline_wgs/scripts/variant_filtering/mask_sv_gts_for_regenotyping.py \
      --quality-field "~{sv_mask_field}" \
    | /opt/pancan_germline_wgs/scripts/variant_filtering/mask_sv_gts_for_regenotyping.py \
      --quality-field "~{sv_mask_field}" \
      --output-vcf sv.for_jaccard.vcf.gz
    for vc in indel sv; do
      tabix -p vcf -f $vc.for_jaccard.vcf.gz
    done
    cat <( echo "#vid" ) union.samples.list | paste -s > jaccard.header
    bcftools concat -a \
      -Oz -o all.for_jaccard.vcf.gz \
      indel.for_jaccard.vcf.gz \
      sv.for_jaccard.vcf.gz
    tabix -p vcf -f all.for_jaccard.vcf.gz
    /opt/pancan_germline_wgs/scripts/utilities/vcf2dosage.py \
      --vcf-in all.for_jaccard.vcf.gz \
    | gzip -c \
    > jaccard.input.tsv.gz
    /opt/pancan_germline_wgs/scripts/variant_filtering/calc_nonref_jaccards.R \
      --genotype-matrix jaccard.input.tsv.gz \
      --pairs-tsv candidate_hits.pairs.tsv.gz \
      --matrix-is-ad \
    | gzip -c \
    > jaccard.output.tsv.gz

    # Step 4: Resolve clusters to assign final variant integrations
    /opt/pancan_germline_wgs/scripts/variant_filtering/filter_indel_sv_clusters.py \
      --pair-distance candidate_hits.pairs.tsv.gz \
      --pair-jaccard jaccard.output.tsv.gz \
      --minimum-jaccard ~{min_nonref_jaccard} \
    | gzip -c \
    > "~{out_fname}"

    # Count number of input indels and SVs in final clusters:
    
    # Total clusters
    zcat "~{out_fname}" | awk '!/^#/' | wc -l > total_clusters.count.txt

    # Number of unique indels
    zcat "~{out_fname}" | awk '!/^#/' \
    | cut -f2 | sed 's/,/\n/g' \
    | sort | uniq | wc -l \
    > indel.variant_count.txt

    # Number of clusters involving multiple indels
    zcat "~{out_fname}" | awk -v FS="\t" '!/^#/ && $2 ~ /,/' | wc -l \
    > multi_indel.clusters.count.txt

    # Number of unique SVs
    zcat "~{out_fname}" | awk '!/^#/' \
    | cut -f1 | sed 's/,/\n/g' \
    | sort | uniq | wc -l \
    > sv.variant_count.txt

    # Number of clusters involving multiple SVs
    zcat "~{out_fname}" | awk -v FS="\t" '!/^#/ && $1 ~ /,/' | wc -l \
    > multi_sv.clusters.count.txt
  >>>

  output {
    File final_clusters = "~{out_fname}"
    File combined_header = "combined.header.vcf.gz"
    Int n_clusters = read_int("total_clusters.count.txt")
    Int n_multi_indel_clusters = read_int("multi_indel.clusters.count.txt")
    Int n_clustered_indels = read_int("indel.variant_count.txt")
    Int n_multi_sv_clusters = read_int("multi_sv.clusters.count.txt")
    Int n_clustered_svs = read_int("sv.variant_count.txt")
  }

  runtime {
    docker: g2c_analysis_docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
    bootDiskSizeGb: 15
    preemptible: 3
    maxRetries: 1
  }
}


# Disambiguate large indels from small SVs and collapse clustered variants
task ResolveClusters {
  input {
    File indel_vcf
    File indel_vcf_idx
    
    File sv_vcf
    File sv_vcf_idx

    File clusters
    File combined_header

    String g2c_analysis_docker

    Float mem_gb = 3.5
    Int n_cpu = 2
  }

  Int sort_mem_mb = floor(1000 * (mem_gb - 2))
  Int disk_gb = ceil(3.5 * size([sv_vcf, indel_vcf], "GB")) + 10
  String output_prefix = basename(sv_vcf, ".svs.vcf.gz")

  command <<<
    set -eu -o pipefail

    # Start heartbeat to avoid silent VM death
    while true; do
      echo "[ResolveClusters] still running at $(date)"
      sleep 60
    done &
    HEARTBEAT_PID=$!
    trap "kill $HEARTBEAT_PID 2>/dev/null || true" EXIT

    # Ensure VCF indexes localize to same directory as VCFs
    if [ "~{indel_vcf}.tbi" != "~{indel_vcf_idx}" ]; then
      cp "~{indel_vcf_idx}" "~{indel_vcf}.tbi"
    fi
    if [ "~{sv_vcf}.tbi" != "~{sv_vcf_idx}" ]; then
      cp "~{sv_vcf_idx}" "~{sv_vcf}.tbi"
    fi
    tabix -p vcf -f "~{combined_header}"

    # Make integrated header for indel & SV VCFs
    tabix -H ~{indel_vcf} | sort -V > indel.header.vcf
    tabix -H ~{sv_vcf} | sort -V > sv.header.vcf
    fgrep "##file" sv.header.vcf > header.vcf
    fgrep "##contig" sv.header.vcf >> header.vcf
    for tag in ALT FILTER INFO FORMAT; do
      cat indel.header.vcf sv.header.vcf | fgrep "##$tag" | sort -V | uniq >> header.vcf
    done
    fgrep "##CPX" sv.header.vcf >> header.vcf
    tabix -H "~{combined_header}" | fgrep -v "##" >> header.vcf

    # Enforce expected sample order for indel and SV VCFs
    bcftools query -l "~{combined_header}" > samples.list
    bcftools view \
      --samples-file samples.list \
      --force-samples \
      -Oz -o input.indels.reordered.vcf.gz \
      ~{indel_vcf}
    tabix -p vcf -f input.indels.reordered.vcf.gz
    rm ~{indel_vcf} 
    bcftools view \
      --samples-file samples.list \
      --force-samples \
      -Oz -o input.svs.reordered.vcf.gz \
      ~{sv_vcf}
    tabix -p vcf -f input.svs.reordered.vcf.gz
    rm ~{sv_vcf} 

    # Integrate indel and SV VCFs
    /opt/pancan_germline_wgs/scripts/variant_filtering/integrate_gatk_vcfs.py \
      --indel-vcf input.indels.reordered.vcf.gz \
      --sv-vcf input.svs.reordered.vcf.gz \
      --clusters "~{clusters}" \
      --out-vcf-header header.vcf \
      --out-prefix "~{output_prefix}.unsorted"
    for vc in indel sv; do
      bcftools sort \
        --max-mem "~{sort_mem_mb}M" \
        -Oz -o "~{output_prefix}.$vc.vcf.gz" \
        "~{output_prefix}.unsorted.$vc.vcf.gz"
      tabix -p vcf -f "~{output_prefix}.$vc.vcf.gz"
      bcftools index -n "~{output_prefix}.$vc.vcf.gz" > $vc.variant_count.txt
    done
  >>>

  output {
    File integrated_indel_vcf = "~{output_prefix}.indel.vcf.gz"
    File integrated_indel_vcf_idx = "~{output_prefix}.indel.vcf.gz.tbi"
    Float integrated_indel_vcf_size = size(integrated_indel_vcf, "GB")
    Int n_integrated_indels = read_int("indel.variant_count.txt")

    File integrated_sv_vcf = "~{output_prefix}.sv.vcf.gz"
    File integrated_sv_vcf_idx = "~{output_prefix}.sv.vcf.gz.tbi"
    Float integrated_sv_vcf_size = size(integrated_sv_vcf, "GB")
    Int n_integrated_svs = read_int("sv.variant_count.txt")
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


# Compute the size of the second-smallest interval in a BED file
# Capped at 100kb for robustness
task DefineLargeSVCutoff {
  input {
    File intervals
    String concat_cmd = "zcat"
    Int bp_buffer = 2
    Int max_cutoff = 100000
    String docker
  }

  Int disk_gb = ceil(2 * size(intervals, "GB")) + 10

  command <<<
    set -eu -o pipefail

    echo 1000000000 > size.txt

    ~{concat_cmd} ~{intervals} \
    | grep -ve '^#' \
    | awk -v FS="\t" '{ print $3-$2 }' \
    | sort -nk1,1 \
    | head -n2 \
    | tail -n1 \
    > size.txt || true
  >>>

  output {
    Int raw_cutoff = read_int("size.txt") - bp_buffer
    Int cutoff = if raw_cutoff > max_cutoff then max_cutoff else raw_cutoff
  }

  runtime {
    docker: docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    maxRetries: 1
  }  
}


# Split large and small SVs from a sorted, integrated VCF
task ExtractLargeSvs {
  input {
    File vcf
    File vcf_idx
    Int size_cutoff
    String bcftools_docker
  }

  Int disk_gb = ceil(3 * size(vcf, "GB")) + 20
  String main_outfile = basename(vcf)
  String large_outfile = basename(vcf, ".vcf.gz") + ".large.vcf.gz"

  command <<<
    set -eu -o pipefail

    # Start heartbeat to avoid silent VM death
    while true; do
      echo "[ExtractLargeSvs] still running at $(date)"
      sleep 60
    done &
    HEARTBEAT_PID=$!
    trap "kill $HEARTBEAT_PID 2>/dev/null || true" EXIT

    echo "Starting bcftools large SV extraction $(date)"
    bcftools view \
      -i 'INFO/SVLEN >= ~{size_cutoff} | INFO/SVTYPE = "CTX"' \
      -Oz -o ~{large_outfile} \
      ~{vcf}
    echo "Finished bcftools large SV extraction $(date)"
    tabix -p vcf -f ~{large_outfile}

    # Check if any large SVs were present in this shard
    # If not, we can simply return the input VCF as the output VCF
    n_large=$( bcftools index -n ~{large_outfile} )
    if [ $n_large -eq 0 ]; then
      rm ~{large_outfile} "~{large_outfile}.tbi"
      mv ~{vcf} ~{main_outfile}
      mv ~{vcf_idx} "~{main_outfile}.tbi"
      exit 0
    fi

    # Reciprocally, we need to *exclude* large SVs from the 
    # input VCF to produce two disjoint output files
    echo "Starting bcftools large SV exclusion $(date)"
    bcftools view \
      -e 'INFO/SVLEN >= ~{size_cutoff} | INFO/SVTYPE = "CTX"' \
      -Oz -o ~{main_outfile} \
      ~{vcf}
    echo "Finished bcftools large SV exclusion $(date)"
    tabix -p vcf -f ~{main_outfile}
  >>>

  output {
    File notlarge_vcf = main_outfile
    File notlarge_vcf_idx = main_outfile + ".tbi"

    File? large_vcf = large_outfile
    File? large_vcf_idx = large_outfile + ".tbi"
  }

  runtime {
    docker: bcftools_docker
    memory: "8 GB"
    cpu: 4
    disks: "local-disk " + disk_gb + " SSD"
    preemptible: 3
    maxRetries: 1
  }  
}


# Small task to write a summary log of workflow inputs and outputs for future reference
task WriteSummaryLog {
  input {
    Int snv_count
    Int small_indel_count
    Int large_indel_count
    Int large_indel_ge50bp_count

    Int eligible_sv_count
    Int passthrough_sv_count
    
    Int cluster_count
    
    Int multi_indel_cluster_count
    Int clustered_indel_count
    
    Int multi_sv_cluster_count
    Int clustered_sv_count
    
    String linux_docker
  }

  Int total_indel_count = small_indel_count + large_indel_count
  Int total_sv_count = eligible_sv_count + passthrough_sv_count

  command <<<
    set -eu -o pipefail

    cat << EOF > "UnifyGatkCallsets.stats.txt"
Input variant counts from GATK-HC
=================================
- SNVs: ~{snv_count}
- All indels: ~{total_indel_count}
- Small indels not eligible for SV clustering: ~{small_indel_count}
- Large indels eligible for SV clustering: ~{large_indel_count}
- Indels 50bp or larger: ~{large_indel_ge50bp_count}

Input variant counts from GATK-SV
=================================
- All SVs: ~{total_sv_count}
- Small CNVs/insertion SVs eligible for indel clustering: ~{eligible_sv_count}
- Other SVs not eligible for indel clustering: ~{passthrough_sv_count}

Indel/SV clustering outcome
===========================
- Final clusters involving indels and SVs: ~{cluster_count}
- Total input indels involved in any cluster: ~{clustered_indel_count}
- Clusters involving two or more indels: ~{multi_indel_cluster_count}
- Total input SVs involved in any cluster: ~{clustered_sv_count}
- Clusters involving two or more SVs: ~{multi_sv_cluster_count}
EOF
  >>>

  output {
    File logfile = "UnifyGatkCallsets.stats.txt"
  }

  runtime {
    docker: linux_docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk 20 HDD"
    preemptible: 3
    maxRetries: 1
  }
}

