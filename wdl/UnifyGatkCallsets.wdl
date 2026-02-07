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

    File genome_file                 # BEDTools-style genome file
    File nmask_bed                   # BED file of N-masked reference intervals
    Int min_interval_size = 100000   # Minimum size of intervals used for SV/indel clustering

    Int min_sv_size = 50
    Float size_scalar = 3            # Maximum fold-difference between sizes of indels and SVs to tolerate

    File snv_partition_intervals
    File indel_partition_intervals
    File sv_partition_intervals

    String g2c_analysis_docker
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
      intervals_bed = MakeClusteringIntervals.intervals_bed,
      disk_gb = ceil(2.2 * EstimateLargeIndelFileSize.sum) + 10,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Partition SVs across clustering intervals
  call ReshardVcfs as SplitSvsForClustering {
    input:
      vcfs = PrepareSvs.eligible_sv_vcf,
      intervals_bed = MakeClusteringIntervals.intervals_bed,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Cluster large indels and SVs
  # TODO: implement this

  # Re-integrate small and large indels
  # TODO: implement this

  # Postprocess SNVs
  call Utils.Sum as EstimateSnvFileSize {
    input:
      values = select_all(SplitGatkHcBySize.snv_vcf_size)
  }
  call ReshardVcfs as PartitionSnvOutputs {
    input:
      vcfs = select_all(SplitGatkHcBySize.snv_vcf),
      intervals_bed = snv_partition_intervals,
      rename = true,
      delete_empty = true,
      disk_gb = ceil(2.2 * EstimateSnvFileSize.sum) + 10,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Postprocess indels
  # TODO: implement this
  # Don't forget: need to add short indels

  # Postprocess SVs
  # TODO: implement this
  # Don't forget: need to add passthrough SVs

  output {
    Array[File?] cleaned_snv_vcfs = PartitionSnvOutputs.sharded_vcfs
    Array[File?] cleaned_snv_vcf_idxs = PartitionSnvOutputs.sharded_vcf_idxs
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
      -i 'INFO/SVTYPE ~ "DEL,DUP,INS" & INFO/SVLEN <= ~{max_size}' \
      -Oz -o "~{elig_outfile}" \
      ~{vcf}
    tabix -p vcf "~{elig_outfile}"

    # Exclusion
    bcftools view \
      -e 'INFO/SVTYPE ~ "DEL,DUP,INS" & INFO/SVLEN <= ~{max_size}' \
      -Oz -o "~{pt_outfile}" \
      ~{vcf}
    tabix -p vcf "~{pt_outfile}"
  >>>

  output {
    File eligible_sv_vcf = "~{elig_outfile}"
    File eligible_sv_vcf_idx = "~{elig_outfile}.tbi"
    File passthrough_sv_vcf = "~{pt_outfile}"
    File passthrough_sv_vcf_idx = "~{pt_outfile}.tbi"
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
    File intervals_bed

    Boolean rename = false
    Boolean delete_empty = false

    String g2c_analysis_docker

    Int disk_gb = 50
    Float mem_gb = 3.5
    Int n_cpu = 2
    Int boot_gb = 25
  }

  String rename_cmd = if rename then "| /opt/pancan_germline_wgs/scripts/qc/vcf_qc/set_g2c_qc_variant_ids.py" else ""

  Int sort_mem_mb = floor(1000 * (mem_gb - 1))

  command <<<
    set -eu -o pipefail

    # Exit with non-zero status if input array is empty
    if [ ~{length(vcfs)} -eq 0 ]; then
      echo "Error: input VCF array is empty. Need to provide at least one VCF. Exiting"
      sleep 60s
      exit 1
    fi

    # Start heartbeat to avoid silent VM death
    (
      while true; do
        echo "[ConcatVcfs] still running at $(date)"
        sleep 60
      done
    ) &
    HEARTBEAT_PID=$!

    # Reshard variants
    cat ~{write_lines(vcfs)} > vcf.inputs.list
    /opt/pancan_germline_wgs/scripts/utilities/reshard_vcfs.py \
      --vcf-list ~{write_lines(vcfs)} \
      --intervals ~{intervals_bed}

    # To reduce disk pressure, we delete the localized copies of raw VCFs
    cat vcf.inputs.list | xargs -I {} rm {}

    # Next, we sort, deduplicate, rename, and reindex each sharded VCF
    mkdir clean_outputs
    while read chrom start end iid; do
      vcf="$iid.vcf.gz"

      # Delete empty VCFs if optioned
      if ~{delete_empty}; then
        if [ $( bcftools index -n $vcf ) -eq 0 ]; then
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
        -D --threads ~{n_cpu} \
        -Oz -o $iid.clean.vcf.gz
      mv $iid.clean.vcf.gz $vcf

      tabix -p vcf -f $vcf

      mv $vcf clean_outputs/
      mv $vcf.tbi clean_outputs/
    done < ~{intervals_bed}

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
    preemptible: 1
    maxRetries: 1
  }
}

