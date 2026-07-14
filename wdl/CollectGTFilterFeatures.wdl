# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Collect features for genotype filtering model for one or more input VCFs
# Given the design of G2C, this workflow assumes input VCFs contain 
# exactly one variant class (SNV, indel, SV)


version 1.0


import "Utilities.wdl" as Utils


workflow CollectGTFilterFeatures {
	input {
    # Two ways to provide VCF information: as arrays for VCF & indexes, or
    # as a two-column .tsv with URIs for VCF and index. If both are provided,
    # the array-style inputs will be used.
    Array[File]? vcf_array
    Array[File]? vcf_idx_array
    File? vcf_info_tsv

    Array[File]? bed_features            # Optional array of BED tracks with features to annotate vs. variant start/end coordinates
    Array[String]? bed_feature_names     # Optional array of feature names for `bed_features` (specified in the same order)
    Array[File]? bigwig_features         # Optional array of bigWig tracks with features to annotate vs. variant start/end coordinates
    Array[String]? bigwig_feature_names  # Optional array of feature names for `bigwig_features` (specified in the same order)

    String output_prefix

    String g2c_analysis_docker
    String linux_docker
  }

  # Determine method of VCF inputs
  if ( !defined(vcf_array) || !defined(vcf_idx_array) ) {
    call Utils.ReadVcfInfo as ExtractVcfInfo {
      input:
        vcf_info = select_first([vcf_info_tsv, '']),
        linux_docker = linux_docker
    }
  }
  Array[File] vcfs = select_first([ExtractVcfInfo.vcf_uris, vcf_array])
  Array[File] vcf_idxs = select_first([ExtractVcfInfo.vcf_tbi_uris, vcf_idx_array])

  # Collect features in parallel for each VCF
  scatter ( vcf_pair in zip(vcfs, vcf_idxs) ) {

    File vcf = vcf_pair.left
    File vcf_idx = vcf_pair.right

    # Collect no-call rates per variant class for each sample
    call CollectNoCallRates {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        g2c_analysis_docker = g2c_analysis_docker
    }

    # Collect site-level features
    call CollectSiteFeatures {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        bed_features = bed_features,
        bed_feature_names = bed_feature_names,
        bigwig_features = bigwig_features,
        bigwig_feature_names = bigwig_feature_names,
        g2c_analysis_docker = g2c_analysis_docker
    }

    # Collect genotype-level features
    call CollectGtFeatures {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }

  # Sum no-call counts
  call SumNoCallCounts {
    input:
      shards = CollectNoCallRates.nocall_counts,
      output_prefix = output_prefix,
      docker = g2c_analysis_docker
  }

  # Concatenate site-level features
  call Utils.ConcatTextFiles as ConcatSiteFeatures {
    input:
      shards = CollectSiteFeatures.site_features,
      concat_command = "zcat",
      compression_command = "gzip",
      input_has_header = true,
      output_filename = output_prefix + ".site_features.tsv.bgz",
      docker = g2c_analysis_docker
  }
  call IndexFeatures as IndexSiteFeatures {
    input:
      features = ConcatSiteFeatures.merged_file,
      docker = g2c_analysis_docker
  }

  # Concatenate genotype-level features
  call Utils.ConcatTextFiles as ConcatGtFeatures {
    input:
      shards = CollectGtFeatures.gt_features,
      concat_command = "zcat",
      compression_command = "gzip",
      input_has_header = true,
      output_filename = output_prefix + ".gt_features.tsv.bgz",
      docker = g2c_analysis_docker
  }
  call IndexFeatures as IndexGtFeatures {
    input:
      features = ConcatGtFeatures.merged_file,
      docker = g2c_analysis_docker
  }

  output {
    File nocall_counts = SumNoCallCounts.nocall_counts

    File site_features = ConcatSiteFeatures.merged_file
    File site_features_idx = IndexSiteFeatures.index

    File gt_features = ConcatGtFeatures.merged_file
    File gt_features_idx = IndexGtFeatures.index
  }
}


# Compute no-call rates for all samples
task CollectNoCallRates {
  input {
    File vcf
    File vcf_idx

    String g2c_analysis_docker
  }

  String outfile = basename(vcf, ".vcf.gz") + ".nocall_counts.tsv.gz"
  Int disk_gb = ceil(1.5 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    # Get list of all samples present in VCF
    bcftools query -l ~{vcf} > samples.list

    # Get total number of records in file
    n_records=$( bcftools index -n ~{vcf} )

    # Tabulate missing genotypes
    bcftools query \
      -i 'INFO/SVTYPE != "MCNV" & GT="mis"' \
      -f '[%SAMPLE\n]' \
      ~{vcf} \
    | cat samples.list - \
    | sort -Vk1,1 \
    | uniq -c \
    | awk -v n="$n_records" -v OFS="\t" '{ print $2, n, $1-1 }' \
    | gzip -c \
    > "~{outfile}"
  >>>

  output {
    File nocall_counts = outfile
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    maxRetries: 1
  }
}


# Collect site-level filtering features
task CollectSiteFeatures {
  input {
    File vcf
    File vcf_idx

    Array[File] bed_features = []
    Array[String] bed_feature_names = []
    Array[File] bigwig_features = []
    Array[String] bigwig_feature_names = []

    String g2c_analysis_docker
  }

  String outfile = basename(vcf, ".vcf.gz") + ".site_features.tsv.bgz"
  Int disk_gb = ceil(2.2 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    # Build options for bed feature annotation
    bfa_cmd=""
    paste \
      ~{write_lines(bed_feature_names)} \
      ~{write_lines(bed_features)} \
    > bed_features.tsv
    while read fname fpath; do
      mv "$fpath" ./
      locpath=$( basename "$fpath" )
      tabix -p bed -f $locpath
      bfa_cmd="$bfa_cmd --feature-bed $fname=$locpath"
    done < bed_features.tsv

    # Build options for bigWig feature annotation
    bwa_cmd=""
    paste \
      ~{write_lines(bigwig_feature_names)} \
      ~{write_lines(bigwig_features)} \
    > bigwig_features.tsv
    while read fname fpath; do
      mv "$fpath" ./
      locpath=$( basename "$fpath" )
      bwa_cmd="$bwa_cmd --feature-bigwig $fname=$locpath"
    done < bigwig_features.tsv

    # Build overall feature collection command
    cmd="/opt/pancan_germline_wgs/scripts/variant_filtering/collect_site_filtering_features.py "
    cmd="$cmd -i \"~{vcf}\" $bfa_cmd $bwa_cmd | bgzip > \"~{outfile}\""
    echo -e "Now collecting site features with the following command:\n$cmd"
    eval "$cmd"

    # Ensure file is complete by generating a tabix index
    tabix -b 2 -e 2 -f -c "#" "~{outfile}"
  >>>

  output {
    File site_features = outfile
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "3.7 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    maxRetries: 1
  }
}


# Collect GT-level filtering features
task CollectGtFeatures {
  input {
    File vcf
    File vcf_idx

    String g2c_analysis_docker
  }

  String outfile = basename(vcf, ".vcf.gz") + ".gt_features.tsv.gz"
  Int disk_gb = ceil(2.2 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    /opt/pancan_germline_wgs/scripts/variant_filtering/collect_gt_filtering_features.py \
      -i "~{vcf}" \
      -o "~{outfile}"
  >>>

  output {
    File gt_features = outfile
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "3.7 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    maxRetries: 1
  }
}


# Sum no-call counts across multiple shards
task SumNoCallCounts {
  input {
    Array[File] shards
    String output_prefix
    String docker
  }

  Int disk_gb = ceil(1.5 * size(shards, "GB")) + 10

  command <<<
    set -eu -o pipefail
    
    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/sum_compressed_distribs.py \
      -o "~{output_prefix}.nocall_counts.tsv" \
      -k 1 \
      ~{sep=" " shards}
    gzip -f "~{output_prefix}.nocall_counts.tsv"
  >>>

  output {
    File nocall_counts = "~{output_prefix}.nocall_counts.tsv.gz"
  }

  runtime {
    docker: docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}


# Make tabix index for feature tsv
task IndexFeatures {
  input {
    File features
    String docker
  }

  Int disk_gb = ceil(1.5 * size(features, "GB")) + 10
  String outfile = basename(features) + ".tbi"

  command <<<
    set -eu -o pipefail

    mv "~{features}" ./

    tabix -b 2 -e 2 -f -c "#" $( basename "~{features}" )
  >>>

  output {
    File index = "~{outfile}"
  }

  runtime {
    docker: docker
    memory: "1.7 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    maxRetries: 1
  }
}
