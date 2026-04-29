# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Nested workflow to reshard an array of VCFs across a predefined set of intervals


version 1.0


import "Utilities.wdl" as Utils


workflow NestedReshardVcfs {
  input {
    # Two ways to provide input VCF information: as arrays for VCF & indexes, or
    # as a two-column .tsv with URIs for VCF and index. If both are provided,
    # the array-style inputs will be used.
    Array[File]? vcfs
    Array[File]? vcf_idxs
    File? vcf_info_tsv

    # Resharding parameters
    File reshard_intervals_bed                     # BED4 file of resharding intervals; fourth column *must* correspond to desired output VCF filename prefix
    Boolean intervals_are_compressed = true
    String? interval_suffix
    File? resharded_vcf_header
    Boolean rename_variants = false
    Boolean keep_empty_resharded_vcfs = false

    # Parallelization parameters
    Boolean shuffle_when_sharding = true           # Should shards be shuffled for random load balancing?
    Int vcfs_per_shard = 10                        # Parallelization control for ReshardVcf tasks
    Int intervals_per_shard = 10                   # Parallelization control for merging resharded VCFs per interval

    Float reshard_task_mem_gb = 7.5
    Int reshard_task_n_cpu = 4

    Float concatenate_task_mem_gb = 3.5
    Int concatenate_task_n_cpu = 2

    String g2c_analysis_docker
    String linux_docker = "ubuntu:plucky-20251001"
  }

  # Determine method of VCF input
  if ( defined(vcfs) && defined(vcf_idxs) ) {
    call Utils.WriteVcfInfo as WriteInputVcfInfo {
      input:
        vcf_uris = select_first([vcfs, []]),
        tbi_uris = select_first([vcf_idxs, []]),
        output_prefix = "input",
        docker = linux_docker
    }
  }
  File input_vcf_info_tsv = select_first(select_all([WriteInputVcfInfo.vcf_info_tsv, vcf_info_tsv]))

  # Chunk input VCFs into N per split
  call Utils.ShardTextFile as ShardInputVcfList {
    input:
      input_file = input_vcf_info_tsv,
      lines_per_split = vcfs_per_shard,
      shuffle = shuffle_when_sharding,
      out_prefix = "input_vcf_chunk",
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Scatter over chunked input VCF list
  scatter ( vcf_info_shard in ShardInputVcfList.shards ) {
    # Read VCF URIs as files
    call Utils.ReadVcfInfo as ReadInputVcfShard {
      input:
        vcf_info = vcf_info_shard,
        linux_docker = linux_docker
    }

    # Run core resharding task
    call Utils.ReshardVcfs {
      input:
        vcfs = ReadInputVcfShard.vcf_uris,
        vcf_idxs = ReadInputVcfShard.vcf_tbi_uris,
        intervals_bed = reshard_intervals_bed,
        intervals_are_compressed = intervals_are_compressed,
        interval_suffix = interval_suffix,
        output_header = resharded_vcf_header,
        rename = rename_variants,
        delete_empty = true,
        mem_gb = reshard_task_mem_gb,
        n_cpu = reshard_task_n_cpu,
        g2c_analysis_docker = g2c_analysis_docker
    }

    # Write tsv with info for resharded VCFs
    Array[File] resharded_chunk_vcfs = select_all(ReshardVcfs.resharded_vcfs)
    Array[File] resharded_chunk_vcf_idxs = select_all(ReshardVcfs.resharded_vcf_idxs)
    if ( length(resharded_chunk_vcfs) > 0 ) {
      call Utils.WriteVcfInfo as WriteReshardedVcfChunk {
        input:
          vcf_uris = resharded_chunk_vcfs,
          tbi_uris = resharded_chunk_vcf_idxs,
          output_prefix = "resharded_vcf_chunk",
          docker = linux_docker
      }
    }
  }

  # Concatenate list of resharded VCFs and group/chunk by intervals
  call ChunkReshardedVcfsByIntervals {
    input:
      vcf_info_tsvs = select_all(WriteReshardedVcfChunk.vcf_info_tsv),
      intervals_bed = reshard_intervals_bed,
      intervals_are_compressed = intervals_are_compressed,
      interval_suffix = interval_suffix,
      intervals_per_chunk = intervals_per_shard,
      shuffle = shuffle_when_sharding,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Scatter over chunks of intervals
  scatter ( interval_vcf_info in ChunkReshardedVcfsByIntervals.vcf_info_chunks ) {
    # Read VCF URIs as files
    call Utils.ReadVcfInfo as ReadIntervalVcfShard {
      input:
        vcf_info = interval_vcf_info,
        linux_docker = linux_docker
    }

    # Localize & concatenate all VCFs per interval
    call ConcatenateIntervalVcfs {
      input:
        vcfs = ReadIntervalVcfShard.vcf_uris,
        vcf_idxs = ReadIntervalVcfShard.vcf_tbi_uris,
        mem_gb = concatenate_task_mem_gb,
        n_cpu = concatenate_task_n_cpu,
        bcftools_docker = g2c_analysis_docker
    }
  }
  Array[File] dense_output_vcfs = select_all(flatten(ConcatenateIntervalVcfs.merged_vcfs))
  Array[File] dense_output_vcf_idxs = select_all(flatten(ConcatenateIntervalVcfs.merged_vcf_idxs))

  # Make empty VCFs for all empty intervals, if any
  if ( ChunkReshardedVcfsByIntervals.n_empty_intervals > 0 && keep_empty_resharded_vcfs ) {
    call MakeEmptyVcfs {
      input:
        empty_interval_names = ChunkReshardedVcfsByIntervals.empty_intervals_list,
        template_vcf = select_first([resharded_vcf_header, select_first(dense_output_vcfs)]),
        bcftools_docker = g2c_analysis_docker
    }
  }

  output {
    Array[File] resharded_vcfs = flatten(select_all([dense_output_vcfs, MakeEmptyVcfs.empty_vcfs]))
    Array[File] resharded_vcf_idxs = flatten(select_all([dense_output_vcf_idxs, MakeEmptyVcfs.empty_vcf_idxs]))
  }
}


task ChunkReshardedVcfsByIntervals {
  input {
    Array[File] vcf_info_tsvs

    File intervals_bed
    Boolean intervals_are_compressed
    String? interval_suffix
    
    Int intervals_per_chunk
    Boolean shuffle

    String g2c_analysis_docker
  }

  String int_cat_cmd = if intervals_are_compressed then "zcat" else "cat"
  String shuffle_cmd = if shuffle then "--shuffle" else ""
  String isuf = if defined(interval_suffix) then ".~{interval_suffix}" else ""

  command <<<
    set -eu -o pipefail

    # Concatenate and sort all resharded VCF URIs
    while read tsv; do 
      cat $tsv
    done < ~{write_lines(vcf_info_tsvs)} \
    | sort -Vk1,1 -k2,2V \
    > input_vcf_info.tsv

    # Determine the list of intervals with at least one input VCF
    cut -f1 input_vcf_info.tsv \
    | xargs -I {} basename {} \
    | sort | uniq \
    > interval_vcf_names.list

    # Determine complement list of intervals with *zero* input VCFs
    ~{int_cat_cmd} ~{intervals_bed} \
    | awk -v isuf="~{isuf}" '{ print $4""isuf".vcf.gz" }' \
    | fgrep -xvf interval_vcf_names.list \
    | sed 's/\.vcf.gz//g' \
    | sort -V | uniq \
    > empty_intervals.list || true

    # Shard intervals
    /opt/pancan_germline_wgs/scripts/utilities/evenSplitter.R \
      -L ~{intervals_per_chunk} \
      ~{shuffle_cmd} \
      interval_vcf_names.list \
      "interval_chunk_"
    find ./ -name "interval_chunk_*" > interval_chunks.list

    # Map VCF URIs to an info tsv for each chunk
    while read member_list; do
      fgrep -f $member_list input_vcf_info.tsv \
      > "vcf_info.$( basename $member_list ).tsv" || true
    done < interval_chunks.list
  >>>

  output {
    Array[File] vcf_info_chunks = glob("vcf_info.*.tsv")
    File empty_intervals_list = "empty_intervals.list"
    Int n_empty_intervals = length(read_lines("empty_intervals.list"))
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk 25 HDD"
    preemptible: 3
    maxRetries: 1
  }
}


# Localize and concatenate all VCFs for each interval
task ConcatenateIntervalVcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs

    String bcftools_docker

    Float mem_gb = 3.5
    Int n_cpu = 2
  }

  Int disk_gb = ceil(3.5 * size(vcfs, "GB")) + 25
  Int sort_mem_mb = floor(1000 * (mem_gb - 2))
  Int concat_threads = floor(2 * n_cpu) - 1

  command <<<
    set -eu -o pipefail

    # Start heartbeat to avoid silent VM death
    (
      while true; do
        echo "[ConcatenateIntervalVcfs] still running at $(date)"
        sleep 60
      done
    ) &
    HEARTBEAT_PID=$!

    # Concatenate and sort all resharded VCF info
    paste \
      ~{write_lines(vcfs)} \
      ~{write_lines(vcf_idxs)} \
    | sort -Vk1,1 -k2,2V \
    > input_vcf_info.tsv

    # Get unique list of intervals to process
    cat input_vcf_info.tsv \
    | cut -f1 \
    | xargs -I {} basename {} \
    | sed 's/\.vcf\.gz//g' \
    | sort -V | uniq \
    > interval_names.list

    # Concatenate all VCFs for each interval
    while read iid; do
      
      # Get list of VCFs to be concatenated
      fgrep "$iid.vcf.gz" input_vcf_info.tsv \
      > $iid.vcf_info.tsv

      # If only one VCF is included, just rename that VCF and continue
      if [ $( wc -l < $iid.vcf_info.tsv ) -eq 1 ]; then
        invcf=$( cut -f1 $iid.vcf_info.tsv | sed -n '1p' )
        intbi=$( cut -f2 $iid.vcf_info.tsv | sed -n '1p' )
        mv "$invcf" $iid.sorted.vcf.gz
        mv "$intbi" $iid.sorted.vcf.gz.tbi
        continue
      fi

      # Otherwise, ensure proper localization of all indexes
      while read vcf tbi; do
        if [ "$tbi" != "$vcf.tbi" ]; then
          cp "$tbi" "$vcf.tbi"
        fi
      done < $iid.vcf_info.tsv

      # Concatenate and sort all records for this interval
      cut -f1 $iid.vcf_info.tsv > $iid.input_vcfs.list
      echo "Processing interval: $iid"
      cat $iid.input_vcfs.list
      bcftools concat -a -D \
        --file-list $iid.input_vcfs.list \
        --threads ~{concat_threads} \
        -Oz -o $iid.concat.vcf.gz
      bcftools index -t $iid.concat.vcf.gz
      echo "  Finished concatenating VCF shards"
      bcftools sort \
        --max-mem "~{sort_mem_mb}M" \
        -Oz -o $iid.sorted.vcf.gz \
        $iid.concat.vcf.gz
      tabix -p vcf -f $iid.sorted.vcf.gz
      echo "  Finished sorting concatenated VCF"
      rm $iid.concat.vcf.gz $iid.concat.vcf.gz.tbi

    done < interval_names.list

    kill $HEARTBEAT_PID
    wait $HEARTBEAT_PID 2>/dev/null || true
  >>>

  output {
    Array[File] merged_vcfs = glob("*.sorted.vcf.gz")
    Array[File] merged_vcf_idxs = glob("*.sorted.vcf.gz.tbi")
  }

  runtime {
    docker: bcftools_docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    maxRetries: 1
  }
}


# Makes header-only VCFs for intervals with no qualifying records
# This ensures the overall workflow output has exactly one VCF 
# for each resharding interval, even if empty
task MakeEmptyVcfs {
  input {
    File empty_interval_names
    File template_vcf
    String bcftools_docker
  }

  command <<<
    set -eu -o pipefail

    bcftools view \
      --header-only \
      -Oz -o header.vcf.gz \
      ~{template_vcf}
    tabix -p vcf -f header.vcf.gz

    mkdir outputs

    while read iid; do
      cp header.vcf.gz "outputs/$iid.sorted.vcf.gz"
      cp header.vcf.gz.tbi "outputs/$iid.sorted.vcf.gz.tbi"
    done < ~{empty_interval_names}
    rm header.vcf.gz ~{template_vcf}
  >>>

  output {
    Array[File] empty_vcfs = glob("outputs/*.sorted.vcf.gz")
    Array[File] empty_vcf_idxs = glob("outputs/*.sorted.vcf.gz.tbi")
  }

  runtime {
    docker: bcftools_docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk 25 HDD"
    preemptible: 3
    maxRetries: 1
  }
}

