# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Unify GATK-HC and GATK-SV callsets while collapsing large indels & small SVs


version 1.0


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

    File genome_file                        # BEDTools-style genome file
    File nmask_bed                          # BED file of N-masked reference intervals
    Int min_interval_size = 100000          # Minimum size of intervals used for SV/indel clustering

    Int min_sv_size = 50
    Float size_scalar = 3                   # Maximum fold-difference between sizes of indels and SVs to tolerate

    # BED4 files of final output intervals for resharding
    # Fourth column must correspond to desired output VCF name / interval name
    File snv_partition_intervals
    File indel_partition_intervals
    File sv_partition_intervals
    Int final_intervals_per_shard = 25      # Parallelization control for final partitioning task
    Float final_partition_disk_scalar = 1.5 # Disk sizing parameter for final partitioning task

    String g2c_analysis_docker
    String g2c_analysis_docker_dev_tmp      # REMOVE THIS ONCE SUCCESSFUL
    String linux_docker = "ubuntu:plucky-20251001"
  }

  Int min_large_indel_size = floor(min_sv_size / size_scalar)

  # Make indel/sv clustering intervals
  call MakeClusteringIntervals {
    input:
      genome_file = genome_file,
      nmask_bed = nmask_bed,
      min_nmask_size = ceil(2 * size_scalar * min_sv_size),
      min_interval_size = min_interval_size,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Determine method of GATK-HC VCF input
  if ( !defined(gatkhc_vcf_array) || !defined(gatkhc_vcf_idx_array) ) {
    call Utils.ExtractVcfArrays as ExtractGatkhcInfo {
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
  call Utils.Max as CalcMaxIndelSize {
    input:
      values = select_all(SplitGatkHcBySize.largest_indel_size)
  }

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

  # Partition large indels across clustering intervals
  call Utils.Sum as EstimateLargeIndelFileSize {
    input:
      values = select_all(SplitGatkHcBySize.large_indel_vcf_size)
  }
  call ReshardVcfs as SplitIndelsForClustering {
    input:
      vcfs = select_all(SplitGatkHcBySize.large_indel_vcf),
      vcf_idxs = select_all(SplitGatkHcBySize.large_indel_vcf_idx),
      intervals_bed = MakeClusteringIntervals.intervals_bed,
      interval_suffix = "large_indels",
      rename = true,
      disk_gb = ceil(2.5 * EstimateLargeIndelFileSize.sum) + 10,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Partition SVs across clustering intervals
  call ReshardVcfs as SplitSvsForClustering {
    input:
      vcfs = PrepareSvs.eligible_sv_vcf,
      vcf_idxs = PrepareSvs.eligible_sv_vcf_idx,
      intervals_bed = MakeClusteringIntervals.intervals_bed,
      interval_suffix = "svs",
      rename = true,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Cluster large indels and SVs
  Int n_cluster_shards = length(SplitIndelsForClustering.sharded_vcfs)
  scatter ( i in range(n_cluster_shards) ) {

    call DefineClusters {
      input:
        indel_vcf = select_first(select_all([SplitIndelsForClustering.sharded_vcfs[i]])),
        indel_vcf_idx = select_first(select_all([SplitIndelsForClustering.sharded_vcf_idxs[i]])),
        sv_vcf = select_first(select_all([SplitSvsForClustering.sharded_vcfs[i]])),
        sv_vcf_idx = select_first(select_all([SplitSvsForClustering.sharded_vcf_idxs[i]])),
        genome_file = genome_file,
        size_scalar = size_scalar,
        g2c_analysis_docker = g2c_analysis_docker
    }

    call ResolveClusters {
      input:
        indel_vcf = select_first(select_all([SplitIndelsForClustering.sharded_vcfs[i]])),
        indel_vcf_idx = select_first(select_all([SplitIndelsForClustering.sharded_vcf_idxs[i]])),
        sv_vcf = select_first(select_all([SplitSvsForClustering.sharded_vcfs[i]])),
        sv_vcf_idx = select_first(select_all([SplitSvsForClustering.sharded_vcf_idxs[i]])),
        clusters = DefineClusters.final_clusters,
        g2c_analysis_docker = g2c_analysis_docker_dev_tmp
    }

    # Possible TODO: could extract records matching clusters here and integrate them in a separate task
    # This would improve speed because the vast majority of indel records could simply be written out

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

  # Postprocess SNVs
  call Utils.Sum as EstimateSnvFileSize {
    input:
      values = select_all(SplitGatkHcBySize.snv_vcf_size)
  }
  call Utils.ShardTextFile as ShardSnvIntervals {
    input:
      input_file = snv_partition_intervals,
      lines_per_split = final_intervals_per_shard,
      out_prefix = "snv_output_partitions",
      g2c_analysis_docker = g2c_analysis_docker
  }
  scatter ( interval_shard in ShardSnvIntervals.shards ) {
    call ReshardVcfs as PartitionSnvOutputs {
      input:
        vcfs = select_all(SplitGatkHcBySize.snv_vcf),
        vcf_idxs = select_all(SplitGatkHcBySize.snv_vcf_idx),
        intervals_bed = interval_shard,
        rename = true,
        delete_empty = true,
        disk_gb = ceil(final_partition_disk_scalar * EstimateSnvFileSize.sum) + 10,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }

  # Postprocess indels
  Array[File] all_indel_vcfs = select_all(flatten([SplitGatkHcBySize.small_indel_vcf,
                                                   ResolveClusters.integrated_indel_vcf]))
  Array[File] all_indel_vcf_idxs = select_all(flatten([SplitGatkHcBySize.small_indel_vcf_idx,
                                                       ResolveClusters.integrated_indel_vcf_idx]))
  call Utils.Sum as EstimateIndelFileSize {
    input:
      values = select_all(flatten([SplitGatkHcBySize.small_indel_vcf_size,
                                   ResolveClusters.integrated_indel_vcf_size]))
  }
  call Utils.ShardTextFile as ShardIndelIntervals {
    input:
      input_file = indel_partition_intervals,
      lines_per_split = final_intervals_per_shard,
      out_prefix = "indel_output_partitions",
      g2c_analysis_docker = g2c_analysis_docker
  }
  scatter ( interval_shard in ShardIndelIntervals.shards ) {
    call ReshardVcfs as PartitionIndelOutputs {
      input:
        vcfs = all_indel_vcfs,
        vcf_idxs = all_indel_vcf_idxs,
        intervals_bed = interval_shard,
        rename = true,
        delete_empty = true,
        disk_gb = ceil(final_partition_disk_scalar * EstimateIndelFileSize.sum) + 10,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }

  # Postprocess SVs
  # Don't forget: need to add passthrough SVs
  Array[File] all_sv_vcfs = select_all(flatten([PrepareSvs.passthrough_sv_vcf,
                                                ResolveClusters.integrated_sv_vcf]))
  Array[File] all_sv_vcf_idxs = select_all(flatten([PrepareSvs.passthrough_sv_vcf_idx,
                                                    ResolveClusters.integrated_sv_vcf_idx]))
  call Utils.Sum as EstimateSvFileSize {
    input:
      values = select_all(flatten([PrepareSvs.passthrough_sv_vcf_size,
                                   ResolveClusters.integrated_sv_vcf_size]))
  }
  call Utils.ShardTextFile as ShardSvIntervals {
    input:
      input_file = sv_partition_intervals,
      lines_per_split = final_intervals_per_shard,
      out_prefix = "sv_output_partitions",
      g2c_analysis_docker = g2c_analysis_docker
  }
  scatter ( interval_shard in ShardSvIntervals.shards ) {
    call ReshardVcfs as PartitionSvOutputs {
      input:
        vcfs = all_sv_vcfs,
        vcf_idxs = all_sv_vcf_idxs,
        intervals_bed = interval_shard,
        rename = true,
        delete_empty = true,
        disk_gb = ceil(final_partition_disk_scalar * EstimateSvFileSize.sum) + 10,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }

  output {
    Array[File?] cleaned_snv_vcfs = select_all(flatten(PartitionSnvOutputs.sharded_vcfs))
    Array[File?] cleaned_snv_vcf_idxs = select_all(flatten(PartitionSnvOutputs.sharded_vcf_idxs))

    Array[File?] cleaned_indel_vcfs = select_all(flatten(PartitionIndelOutputs.sharded_vcfs))
    Array[File?] cleaned_indel_vcf_idxs = select_all(flatten(PartitionIndelOutputs.sharded_vcf_idxs))

    Array[File?] cleaned_sv_vcfs = select_all(flatten(PartitionSvOutputs.sharded_vcfs))
    Array[File?] cleaned_sv_vcf_idxs = select_all(flatten(PartitionSvOutputs.sharded_vcf_idxs))

    File cluster_assignment_log = ConcatClusterLogs.merged_file
  }
}


# Make intervals for indel/SV clustering based on N-masked reference regions
task MakeClusteringIntervals {
  input {
    File genome_file
    File nmask_bed
    Int min_nmask_size
    Int min_interval_size
    String g2c_analysis_docker
  }

  command <<<
    set -eu -o pipefail

    # Prep N-mask
    zcat ~{nmask_bed} \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | bedtools merge -i - \
    | awk -v ms="~{min_nmask_size}" -v FS="\t" -v OFS="\t" \
      '{ if ($3-$2>=ms) print }' \
    | awk -v md="~{min_interval_size}" -v FS="\t" -v OFS="\t" \
      'NR==1 { prev_chr=$1; prev_end=$3; print; next }
      { if ($1 != prev_chr || $2 >= prev_end + md) { prev_chr=$1; prev_end=$3; print } }' \
    | bgzip -c \
    > nmask.filtered.bed.gz

    # Subtract qualifying N-masked intervals from genome file
    awk -v FS="\t" -v OFS="\t" '{ print $1, 0, $2 }' ~{genome_file} \
    | bedtools subtract -a - -b nmask.filtered.bed.gz \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | awk -v OFS="\t" '{ print $0, "clustering_interval_"NR }' \
    | bgzip -c \
    > clustering.intervals.bed.gz
  >>>

  output {
    File intervals_bed = "clustering.intervals.bed.gz"
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

  Int disk_gb = ceil(3 * size(vcf, "GB")) + 10
  String out_prefix = basename(vcf, ".vcf.gz")
  String snv_out_vcf = out_prefix + ".snv.vcf.gz"
  String si_out_vcf = out_prefix + ".small_indel.vcf.gz"
  String li_out_vcf = out_prefix + ".large_indel.vcf.gz"

  command <<<
    set -eu -o pipefail

    # Ensure tabix index is localized to the same directory as the input VCF
    if [ "~{vcf}.tbi" != "~{vcf_idx}" ]; then
      mv "~{vcf_idx}" "~{vcf}.tbi"
    fi 

    # Split input VCF
    echo "0" > "~{out_prefix}.largest_indel_size.txt"
    /opt/pancan_germline_wgs/scripts/gatkhc_helpers/split_gatkhc_by_size.py \
      -i ~{vcf} \
      --large-indel-min-size ~{large_indel_size} \
      --largest-indel-log "~{out_prefix}.largest_indel_size.txt" \
      -o ~{out_prefix}

    # Tabix non-empty output VCFs and delete empty VCFs
    while read subvcf; do
      tabix -p vcf -f $subvcf
      if [ $( bcftools index -n $subvcf ) -eq 0 ]; then
        rm $subvcf $subvcf.tbi
      fi
    done < <( find ./ -name "~{out_prefix}.*.vcf.gz" )

    # Validation check: since all outputs are optional, this task can seemingly
    # be interpreted by Cromwell as a success even if it errors out. Thus,
    # we need to explicitly exit non-zero if there are no VCFs in pwd
    if [ ~{strict} ]; then
      if [ $( find ./ -name "~{out_prefix}.*.vcf.gz" | wc -l ) -lt 1 ]; then
        echo "Likely error: no output VCFs were generated. This is usually unexpected and indicates a silent error."
        sleep 60s
        exit 1
      fi
    fi
  >>>

  output {
    File? snv_vcf = snv_out_vcf
    File? snv_vcf_idx = snv_out_vcf + ".tbi"
    Float? snv_vcf_size = size(snv_out_vcf, "GB")

    File? small_indel_vcf = si_out_vcf
    File? small_indel_vcf_idx = si_out_vcf + ".tbi"
    Float? small_indel_vcf_size = size(si_out_vcf, "GB")

    File? large_indel_vcf = li_out_vcf
    File? large_indel_vcf_idx = li_out_vcf + ".tbi"
    Float? large_indel_vcf_size = size(li_out_vcf, "GB")

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

    # Exclusion
    bcftools view \
      -e 'INFO/SVTYPE = "DEL,DUP,INS" & INFO/SVLEN <= ~{max_size}' \
      -Oz -o "~{pt_outfile}" \
      ~{vcf}
    tabix -p vcf "~{pt_outfile}"
  >>>

  output {
    File eligible_sv_vcf = "~{elig_outfile}"
    File eligible_sv_vcf_idx = "~{elig_outfile}.tbi"

    File passthrough_sv_vcf = "~{pt_outfile}"
    File passthrough_sv_vcf_idx = "~{pt_outfile}.tbi"
    Float passthrough_sv_vcf_size = size(passthrough_sv_vcf, "GB")
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


# Reshard one or more input VCFs across prespecified intervals
task ReshardVcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    File intervals_bed
    Boolean intervals_are_compressed = true
    String? interval_suffix

    Boolean rename = false
    Boolean delete_empty = false

    String g2c_analysis_docker

    Int disk_gb = 50
    Float mem_gb = 7.5
    Int n_cpu = 4
    Int boot_gb = 25
    Int n_preemptible = 1
    Int ulimit = 4096
  }

  String rename_cmd = if rename then "| /opt/pancan_germline_wgs/scripts/qc/vcf_qc/set_g2c_qc_variant_ids.py" else ""
  String int_cat_cmd = if intervals_are_compressed then "zcat" else "cat"
  String int_bgzip_cmd = if intervals_are_compressed then "| bgzip -c" else ""
  String int_bed_loc = basename(intervals_bed)

  Int sort_mem_mb = floor(1000 * (mem_gb - 2))

  command <<<
    set -eu -o pipefail

    # Update and report ulimit for debugging purposes
    ulimit -n ~{ulimit} || true
    echo -e "\nVM ulimit: $( ulimit -n )\n"

    # Exit with non-zero status if input array is empty
    if [ ~{length(vcfs)} -eq 0 ]; then
      echo "Error: input VCF array is empty. Need to provide at least one VCF. Exiting"
      sleep 60s
      exit 1
    fi

    # Start heartbeat to avoid silent VM death
    (
      while true; do
        echo "[ReshardVcfs] still running at $(date)"
        sleep 60
      done
    ) &
    HEARTBEAT_PID=$!

    # Relocate intervals to pwd, adding a suffix if optioned
    if ~{defined(interval_suffix)}; then
      ~{int_cat_cmd} ~{intervals_bed} \
      | awk -v suf="~{interval_suffix}" -v FS="\t" -v OFS="\t" \
        '{ print $1, $2, $3, $4"."suf }' \
      ~{int_bgzip_cmd} \
      > ~{int_bed_loc}
    else
      cp ~{intervals_bed} ~{int_bed_loc}
    fi

    # Relocate VCF indexes to ensure they match their corresponding VCF
    cat ~{write_lines(vcfs)} > vcf.inputs.list
    cat ~{write_lines(vcf_idxs)} > tbi.inputs.list
    while read vcf; do
      exp_tbi="$vcf.tbi"
      if [ $( fgrep -w $exp_tbi tbi.inputs.list | wc -l ) -eq 0 ]; then
        tbi_base="$( basename $vcf ).tbi"
        cp $( fgrep $tbi_base $tbi.inputs.list ) $exp_tbi
      fi
    done < vcf.inputs.list

    # Reshard variants
    /opt/pancan_germline_wgs/scripts/utilities/reshard_vcfs.py \
      --vcf-list vcf.inputs.list \
      --intervals ~{int_bed_loc}

    # To reduce disk pressure, we delete the localized copies of raw VCFs
    xargs -a vcf.inputs.list rm

    # Next, we sort, deduplicate, rename, and reindex each sharded VCF
    mkdir clean_outputs
    ~{int_cat_cmd} ~{int_bed_loc} > intervals.tmp
    while read chrom start end iid; do
      vcf="$iid.vcf.gz"

      # Delete empty VCFs if optioned
      if ~{delete_empty}; then
        if [ $( bcftools query -f '%ID\n' $vcf | sed -n '10p' | wc -l ) -eq 0 ]; then
          echo -e "$vcf contains no records; removing because `delete_empty` is `true`..."
          rm $vcf
          continue
        fi
      fi

      # Sort, deduplicate, and rename records
      bcftools sort \
        --max-mem "~{sort_mem_mb}M" \
        $vcf \
      ~{rename_cmd} \
      | bcftools norm \
        -d exact \
        --threads ~{n_cpu} \
        -Oz -o $iid.clean.vcf.gz
      mv $iid.clean.vcf.gz $vcf

      tabix -p vcf -f $vcf

      mv $vcf clean_outputs/
      mv $vcf.tbi clean_outputs/
    done < intervals.tmp

    kill $HEARTBEAT_PID
    wait $HEARTBEAT_PID 2>/dev/null || true
  >>>

  output {
    Array[File?] sharded_vcfs = glob("clean_outputs/*vcf.gz")
    Array[File?] sharded_vcf_idxs = glob("clean_outputs/*vcf.gz.tbi")
  }

  runtime {
    docker: g2c_analysis_docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
    bootDiskSizeGb: boot_gb
    preemptible: n_preemptible
    maxRetries: 1
  }
}


# Define clusters of indels and SVs to be integrated
task DefineClusters {
  input {
    File indel_vcf
    File indel_vcf_idx
    
    File sv_vcf
    File sv_vcf_idx

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

    # Ensure VCF indexes localize to same directory as VCFs
    if [ "~{indel_vcf}.tbi" != "~{indel_vcf_idx}" ]; then
      cp "~{indel_vcf_idx}" "~{indel_vcf}.tbi"
    fi
    if [ "~{sv_vcf}.tbi" != "~{sv_vcf_idx}" ]; then
      cp "~{sv_vcf_idx}" "~{sv_vcf}.tbi"
    fi

    # Step 0: if zero indels or SVs are present, there's nothing to do
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
      if [ $( cat $vc.partitioned.sites.list | wc -l ) -gt 1 ]; then
        cat \
          <( zcat $( sed -n '1p' $vc.partitioned.sites.list ) | sed -n '1p' ) \
          <( cat $vc.partitioned.sites.list | xargs -I {} zcat {} \
             | grep -ve '^#' | sort -Vk1,1 -k2,2n -k3,3n -k4,4V ) \
        | bgzip -c \
        > input.$vc.sites.bed.gz
      else
        cp $( sed -n '1p' $vc.partitioned.sites.list ) input.$vc.sites.bed.gz
      fi
      tabix -p bed -f input.$vc.sites.bed.gz
    done

    # Step 2: find candidate overlapping indels & SVs
    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/compare_sites.py \
      -a input.sv.sites.bed.gz \
      -b input.indel.sites.bed.gz \
      -g ~{genome_file} \
      -o candidate_hits \
      --mode both \
      --max-overlap-size-diff ~{size_scalar} \
      --overlap-pad 1 \
      --one-to-many \
      --no-reverse \
      --gzip
    zcat candidate_hits.loj.sites.bed.gz \
    | grep -ve '^#' \
    | awk -v FS="\t" -v OFS="\t" '{ if ($9!="NA") print "sv_"$4, "indel_"$9, $NF }' \
    | sort -Vk1,1 -k3,3n -k2,2V \
    | uniq \
    | gzip -c \
    > candidate_hits.pairs.tsv.gz
    zcat candidate_hits.pairs.tsv.gz \
    | awk -v OFS="\n" '{ print $1, $2 }' \
    | sort -V \
    | uniq \
    > candidate_hits.vids.list

    # Step 3: Compute non-ref GT Jaccard indexes for all candidate pairs
    bcftools query -l ~{indel_vcf} > indel.samples.list
    bcftools query -l ~{sv_vcf} > sv.samples.list
    cat indel.samples.list sv.samples.list | sort -V | uniq > union.samples.list
    for vc in indel sv; do
      case "$vc" in
        "indel")
          invcf="~{indel_vcf}"
          ;;
        "sv")
          invcf="~{sv_vcf}"
          ;;
      esac
      bcftools view \
        --samples-file union.samples.list \
        $invcf \
      | bcftools annotate \
        --set-id "$vc""_%ID" \
      | bcftools view \
        --include 'ID=@candidate_hits.vids.list' \
        -Oz -o $vc.for_jaccard.vcf.gz
      tabix -p vcf -f $vc.for_jaccard.vcf.gz
    done
    cat \
      <( echo "#vid" ) \
      union.samples.list \
    | paste -s \
    > jaccard.header
    bcftools concat -a \
      -Oz -o all.for_jaccard.vcf.gz \
      indel.for_jaccard.vcf.gz \
      sv.for_jaccard.vcf.gz
    tabix -p vcf -f all.for_jaccard.vcf.gz
    bcftools query \
      -f '%ID[\t%GT]\n' \
      all.for_jaccard.vcf.gz \
    | cat jaccard.header - \
    | gzip -c \
    > jaccard.input.tsv.gz
    /opt/pancan_germline_wgs/scripts/variant_filtering/calc_nonref_jaccards.R \
      --genotype-matrix jaccard.input.tsv.gz \
      --pairs-tsv candidate_hits.pairs.tsv.gz \
    | gzip -c \
    > jaccard.output.tsv.gz

    # Step 4: Resolve clusters to assign final variant integrations
    /opt/pancan_germline_wgs/scripts/variant_filtering/filter_indel_sv_clusters.py \
      --pair-distance candidate_hits.pairs.tsv.gz \
      --pair-jaccard jaccard.output.tsv.gz \
      --minimum-jaccard ~{min_nonref_jaccard} \
    | gzip -c \
    > "~{out_fname}"
  >>>

  output {
    File final_clusters = "~{out_fname}"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
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

    String g2c_analysis_docker

    Float mem_gb = 7.5
    Int n_cpu = 4
  }

  Int sort_mem_mb = floor(1000 * (mem_gb - 2))
  Int disk_gb = ceil(2.5 * size([sv_vcf, indel_vcf], "GB")) + 20
  String output_prefix = basename(sv_vcf, ".svs.vcf.gz")

  command <<<
    set -eu -o pipefail

    # Ensure VCF indexes localize to same directory as VCFs
    if [ "~{indel_vcf}.tbi" != "~{indel_vcf_idx}" ]; then
      cp "~{indel_vcf_idx}" "~{indel_vcf}.tbi"
    fi
    if [ "~{sv_vcf}.tbi" != "~{sv_vcf_idx}" ]; then
      cp "~{sv_vcf_idx}" "~{sv_vcf}.tbi"
    fi

    # Make integrated header for indel & SV VCFs
    tabix -H ~{indel_vcf} | sort -V > indel.header.vcf
    tabix -H ~{sv_vcf} | sort -V > sv.header.vcf
    fgrep "##file" sv.header.vcf > header.vcf
    fgrep "##contig" sv.header.vcf >> header.vcf
    for tag in ALT FILTER INFO FORMAT; do
      cat indel.header.vcf sv.header.vcf | fgrep "##$tag" | sort -V | uniq >> header.vcf
    done
    fgrep "##CPX" sv.header.vcf >> header.vcf
    tabix -H all.for_jaccard.vcf.gz | fgrep -v "##" >> header.vcf

    # Step 5: Integrate indel and SV VCFs
    /opt/pancan_germline_wgs/scripts/variant_filtering/integrate_gatk_vcfs.py \
      --indel-vcf ~{indel_vcf} \
      --sv-vcf ~{sv_vcf} \
      --clusters "~{clusters}" \
      --out-vcf-header header.vcf \
      --out-prefix "~{output_prefix}.unsorted"
    for vc in indel sv; do
      bcftools sort \
        --max-mem "~{sort_mem_mb}M" \
        -Oz -o "~{output_prefix}.$vc.vcf.gz" \
        "~{output_prefix}.unsorted.$vc.vcf.gz"
      tabix -p vcf -f "~{output_prefix}.$vc.vcf.gz"
    done
  >>>

  output {
    File integrated_indel_vcf = "~{output_prefix}.indel.vcf.gz"
    File integrated_indel_vcf_idx = "~{output_prefix}.indel.vcf.gz.tbi"
    Float integrated_indel_vcf_size = size(integrated_indel_vcf, "GB")

    File integrated_sv_vcf = "~{output_prefix}.sv.vcf.gz"
    File integrated_sv_vcf_idx = "~{output_prefix}.sv.vcf.gz.tbi"
    Float integrated_sv_vcf_size = size(integrated_sv_vcf, "GB")
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
