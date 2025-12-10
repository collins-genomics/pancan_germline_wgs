# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Post hoc normalization & variant counting from a GATK-HC joint genotyped VCF

# Expected input follows the conventions of GnarlyJointGenotypingPart1.wdl


version 1.0


import "Utilities.wdl" as Utils


workflow PosthocCleanupPart1 {
  input {
    # Input can be specified in two ways: 
    # 1. Arrays of VCFs & tabix indexes:
    Array[File]? vcfs
    Array[File]? vcf_idxs
    # 2. Or a two-column .tsv of VCF and index GCS URIs
    File? vcf_info_tsv

    # Parallelization control (makes workflow metadata easier to query)
    Int vcfs_per_shard = 10

    File samples_list
    String output_prefix

    File ref_fasta
    
    String linux_docker
    String bcftools_docker
    String g2c_analysis_docker
  }

  # If inputs are defined as arrays, write as .tsv for subsequent chunking
  if ( ! defined(vcf_info_tsv) ) {
    call Utils.VcfArrayToTsv {
      input:
        vcfs = select_first([vcfs, []]),
        vcf_idxs = select_first([vcf_idxs, []]),
        output_prefix = output_prefix,
        docker = linux_docker
    }
  }

  # Shard VCF URI list
  call Utils.ShardTextFile as ShardVcfList {
    input:
      input_file = select_first([vcf_info_tsv, VcfArrayToTsv.vcf_info_tsv]),
      lines_per_split = vcfs_per_shard,
      out_prefix = output_prefix + ".",
      g2c_analysis_docker = g2c_analysis_docker
  }
  Int n_chunks = length(ShardVcfList.shards)

  # Scatter over sharded VCF lists
  scatter ( i in range(n_chunks) ) {
    # Extract URIs from sharded file as array of strings and indexes
    call Utils.ExtractVcfArrays {
      input:
        vcf_info = ShardVcfList.shards[i],
        linux_docker = g2c_analysis_docker
    }

    Array[Pair[File, File]] inner_vcf_infos = zip(ExtractVcfArrays.vcf_uris, 
                                                  ExtractVcfArrays.vcf_tbi_uris)
    scatter ( j in range(length(inner_vcf_infos)) ) {

      String out_vcf_fname = output_prefix + "." + i + "." + j + ".norm.vcf.gz"
      File in_vcf_shard = inner_vcf_infos[i].left
      File in_vcf_shard_idx = inner_vcf_infos[i].right

      # Additional cleanup step added for G2C to coerce to parsimonious format
      call NormalizeShortVariants as NormalizeVcf {
        input:
          vcf = in_vcf_shard,
          vcf_idx = in_vcf_shard_idx,
          samples_list = samples_list,
          outfile_name = out_vcf_fname,
          ref_fasta = ref_fasta,
          bcftools_docker = bcftools_docker
      }

      # Count number of non-reference SNVs, insertions, and deletions per sample
      call CountShortVariantsPerSample {
        input:
          vcf = NormalizeVcf.norm_vcf,
          vcf_idx = NormalizeVcf.norm_vcf_idx,
          bcftools_docker = bcftools_docker
      }
    }

    # Inner layer of sum variant counts
    call Utils.SumSvCountsPerSample as SumCountsInner {
      input:
        count_tsvs = CountShortVariantsPerSample.counts_tsv,
        output_prefix = output_prefix + "." + i + ".norm",
        docker = g2c_analysis_docker
    }
  }

  # Outer layer of sum variant counts
  call Utils.SumSvCountsPerSample as SumCounts {
    input:
      count_tsvs = SumCountsInner.summed_tsv,
      output_prefix = output_prefix + ".norm",
      docker = g2c_analysis_docker
  }

  output {
    Array[File] normalized_vcfs = flatten(NormalizeVcf.norm_vcf)
    Array[File] normalized_vcf_idxs = flatten(NormalizeVcf.norm_vcf_idx)
    File counts_per_sample = SumCounts.summed_tsv
  }
}


task NormalizeShortVariants {
  input {
    File vcf
    File vcf_idx
    String outfile_name
    File samples_list
    File ref_fasta
    Float mem_gb = 7.5
    String bcftools_docker
  }

  Int disk_gb = ceil(3 * size([vcf, ref_fasta], "GB")) + 10

  command <<<
    set -eu -o pipefail

    bcftools view \
      --samples-file ~{samples_list} \
      --force-samples \
      ~{vcf} \
    | bcftools annotate \
      -x FORMAT/PGT,FORMAT/PID,FORMAT/PS,INFO/DB \
      -i 'FORMAT/DP>10 & FORMAT/GQ>20 & GT="alt"' \
    | bcftools norm \
      --fasta-ref ~{ref_fasta} \
      --check-ref s \
      --multiallelics - \
      --threads 4 \
      --site-win 100 \
    | bcftools view \
      --include 'INFO/AC>0 & FORMAT/DP>10 & FORMAT/GQ>20 & GT="alt" & alt[0] != "*"' \
      -Oz -o ~{outfile_name}

    tabix -p vcf ~{outfile_name}
  >>>

  output {
    File norm_vcf = outfile_name
    File norm_vcf_idx = "~{outfile_name}.tbi"
  }

  runtime {
    docker: bcftools_docker
    memory: "~{mem_gb} GB"
    cpu: 4
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    maxRetries: 1
  }
}


task CountShortVariantsPerSample {
  input {
    File vcf
    File vcf_idx
    String bcftools_docker
  }

  String outfile = basename(vcf, ".vcf.gz") + ".variant_counts.tsv"
  Int disk_gb = ceil(1.2 * size([vcf], "GB")) + 10

  command <<<
    set -eu -o pipefail

    echo -e "sample\ttype\tcount" > ~{outfile}

    bcftools query \
      --include 'GT == "alt"' \
      --format '[%SAMPLE\t%REF\t%ALT\n]' \
      ~{vcf} \
    | awk -v FS="\t" -v OFS="\t" \
      '{ vl=length($3)-length($2); \
         if (vl>0) { vt="INS" }else \
         if (vl<0) { vt="DEL" }else
         { vt="SNV" };
         print $1, vt }' \
    | sort -Vk1,1 -k2,2V | uniq -c \
    | awk -v OFS="\t" '{ print $2, $3, $1 }' \
    >> ~{outfile}
  >>>

  output {
    File counts_tsv = outfile
  }

  runtime {
    docker: bcftools_docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    maxRetries: 1
  }
}
