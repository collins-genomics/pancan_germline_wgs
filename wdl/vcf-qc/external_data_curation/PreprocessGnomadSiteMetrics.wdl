# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Collect variant summary metrics for gnomAD to match G2C QC
# Designed to run on a single chromosome (based on how gnomAD data are stored online)


version 1.0


import "https://raw.githubusercontent.com/collins-genomics/pancan_germline_wgs/refs/heads/posthoc_qc/wdl/vcf-qc/QcTasks.wdl" as QcTasks


workflow PreprocessGnomadSiteMetrics {
  input {
    File? snv_vcf
    File? snv_vcf_idx
    File? snv_scatter_intervals
    Int? snv_n_samples

    File? sv_vcf
    File? sv_vcf_idx
    File? sv_scatter_intervals
    Int? sv_n_samples

    Float? min_af_bin

    String output_prefix

    String linux_docker
    String g2c_analysis_docker
    String g2c_pipeline_docker
  }

  # Handle optional inputs
  Boolean do_snv = defined(snv_vcf) &&
                   defined(snv_vcf_idx) &&
                   defined(snv_scatter_intervals) &&
                   defined(snv_n_samples)
  Boolean do_sv = defined(sv_vcf) &&
                  defined(sv_vcf_idx) &&
                  defined(sv_scatter_intervals) &&
                  defined(sv_n_samples)

  # Reconcile differences in minimum possible AF between SNVs & SVs, if necessary
  if (do_snv || do_sv) {
    call SetDefaultMinAf {
      input:
        snv_n_samples = snv_n_samples,
        sv_n_samples = sv_n_samples,
        docker = g2c_analysis_docker
    }
  }
  Float min_af_use = select_first([min_af_bin, SetDefaultMinAf.default_min_af])

  # Process SNVs, if provided
  if ( do_snv ) {
    File use_snv_vcf = select_first([snv_vcf])
    File use_snv_vcf_idx = select_first([snv_vcf_idx])
    File use_snv_scatter_intervals = select_first([snv_scatter_intervals])
    Int use_snv_n_samples = select_first([snv_n_samples])

    # Clean SNV scatter intervals
    call QcTasks.ParseIntervals as MakeSnvIntervals {
      input:
        intervals_list = use_snv_scatter_intervals,
        docker = linux_docker
    }
    Array[Array[String]] snv_interval_infos = MakeSnvIntervals.interval_info

    # Scatter SNV VCF over intervals
    scatter ( snv_interval_info in snv_interval_infos ) {

      String snv_interval_coords = snv_interval_info[1]
      String snv_shard_prefix = output_prefix + ".snv." + snv_interval_info[0]

      # Slice gnomAD SNV VCF to interval and drop unnecessary annotations
      call SliceVcf as SliceSnvShard {
        input:
          vcf = use_snv_vcf,
          vcf_idx = use_snv_vcf_idx,
          interval = snv_interval_coords,
          out_prefix = snv_shard_prefix,
          g2c_pipeline_docker = g2c_pipeline_docker
      }

      # Collect site metrics for SNVs
      call QcTasks.CollectSiteMetrics as CollectSnvMetrics {
        input:
          vcf = SliceSnvShard.cleaned_vcf,
          vcf_idx = SliceSnvShard.cleaned_vcf_idx,
          n_samples = use_snv_n_samples,
          min_af_bin = min_af_use,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }

    # Concatenate sites BEDs for SNVs
    call QcTasks.ConcatTextFiles as ConcatSnvs {
      input:
        shards = select_all(CollectSnvMetrics.snv_sites),
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = output_prefix + ".snv.sites.bed.gz",
        docker = g2c_pipeline_docker
    }

    # Collapse indel metrics from SNV input (if any)
    Array[File] short_indel_beds = select_all(CollectSnvMetrics.indel_sites)
    if ( length(short_indel_beds) > 0 ) {
      call QcTasks.ConcatTextFiles as ConcatShortIndels {
        input:
          shards = short_indel_beds,
          concat_command = "zcat",
          sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
          compression_command = "bgzip -c",
          input_has_header = true,
          output_filename = output_prefix + ".short_indel.sites.bed.gz",
          docker = g2c_pipeline_docker
      }
    }

    # Collapse SV metrics from SNV input (if any)
    Array[File] short_sv_beds = select_all(CollectSnvMetrics.sv_sites)
    if ( length(short_sv_beds) > 0 ) {
      call QcTasks.ConcatTextFiles as ConcatShortSvs {
        input:
          shards = short_sv_beds,
          concat_command = "zcat",
          sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
          compression_command = "bgzip -c",
          input_has_header = true,
          output_filename = output_prefix + ".short_sv.sites.bed.gz",
          docker = g2c_pipeline_docker
      }
    }
  }

  # Process SVs, if provided
  if ( do_sv ) {
    File use_sv_vcf = select_first([sv_vcf])
    File use_sv_vcf_idx = select_first([sv_vcf_idx])
    File use_sv_scatter_intervals = select_first([sv_scatter_intervals])
    Int use_sv_n_samples = select_first([sv_n_samples])

    # Clean SV scatter intervals
    call QcTasks.ParseIntervals as MakeSvIntervals {
      input:
        intervals_list = use_sv_scatter_intervals,
        docker = linux_docker
    }
    Array[Array[String]] sv_interval_infos = MakeSvIntervals.interval_info

    # Scatter SV VCF over intervals
    scatter ( sv_interval_info in sv_interval_infos ) {

      String sv_interval_coords = sv_interval_info[1]
      String sv_shard_prefix = output_prefix + ".sv." + sv_interval_info[0]

      # Slice gnomAD SV VCF to interval and drop unnecessary annotations
      call SliceVcf as SliceSvShard {
        input:
          vcf = use_sv_vcf,
          vcf_idx = use_sv_vcf_idx,
          interval = sv_interval_coords,
          out_prefix = sv_shard_prefix,
          g2c_pipeline_docker = g2c_pipeline_docker
      }

      # Collect site metrics for SVs
      call QcTasks.CollectSiteMetrics as CollectSvMetrics {
        input:
          vcf = SliceSvShard.cleaned_vcf,
          vcf_idx = SliceSvShard.cleaned_vcf_idx,
          n_samples = use_sv_n_samples,
          min_af_bin = min_af_use,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }

    # Collapse indel metrics from SV input (if any)
    Array[File] long_indel_beds = select_all(CollectSvMetrics.indel_sites)
    if ( length(long_indel_beds) > 0 ) {
      call QcTasks.ConcatTextFiles as ConcatLongIndels {
        input:
          shards = long_indel_beds,
          concat_command = "zcat",
          sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
          compression_command = "bgzip -c",
          input_has_header = true,
          output_filename = output_prefix + ".long_indel.sites.bed.gz",
          docker = g2c_pipeline_docker
      }
    }

    # Collapse SV metrics from SV input
    call QcTasks.ConcatTextFiles as ConcatLongSvs {
      input:
        shards = select_all(CollectSvMetrics.sv_sites),
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = output_prefix + ".long_sv.sites.bed.gz",
        docker = g2c_pipeline_docker
    }
  }

  # Concatenate sites BEDs for indels
  Array[File] indel_beds = select_all([ConcatShortIndels.merged_file, ConcatLongIndels.merged_file])
  if ( length(indel_beds) > 0 ) {
    call QcTasks.ConcatTextFiles as ConcatIndels {
      input:
        shards = indel_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = output_prefix + ".indel.sites.bed.gz",
        docker = g2c_pipeline_docker
    }    
  }

  # Concatenate sites BEDs for SVs
  Array[File] sv_beds = select_all([ConcatShortSvs.merged_file, ConcatLongSvs.merged_file])
  if ( length(sv_beds) > 0 ) {
    call QcTasks.ConcatTextFiles as ConcatSvs {
      input:
        shards = sv_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = output_prefix + ".sv.sites.bed.gz",
        docker = g2c_pipeline_docker
    }
  }

  # Collapse size distributions
  Array[File] short_size_distribs = select_first([CollectSnvMetrics.size_distrib, []])
  Array[File] sv_size_distribs = select_first([CollectSvMetrics.size_distrib, []])
  Array[File] size_distrib_shards = flatten([short_size_distribs, sv_size_distribs])
  if ( length(size_distrib_shards) > 0 ) {
    call QcTasks.SumCompressedDistribs as SumSizeDistribs {
      input:
        distrib_tsvs = size_distrib_shards,
        out_prefix = output_prefix + ".size_distribution",
        g2c_analysis_docker = g2c_analysis_docker
    }    
  }

  # Collapse AF distributions
  Array[File] short_af_distribs = select_first([CollectSnvMetrics.af_distrib, []])
  Array[File] sv_af_distribs = select_first([CollectSvMetrics.af_distrib, []])
  Array[File] af_distrib_shards = flatten([short_af_distribs, sv_af_distribs])
  if ( length(size_distrib_shards) > 0 ) {
    call QcTasks.SumCompressedDistribs as SumAfDistribs {
      input:
        distrib_tsvs = af_distrib_shards,
        out_prefix = output_prefix + ".af_distribution",
        g2c_analysis_docker = g2c_analysis_docker
    }    
  }

  # Collapse 2D size vs. AF distributions
  Array[File] short_joint_distribs = select_first([CollectSnvMetrics.size_vs_af_distrib, []])
  Array[File] sv_joint_distribs = select_first([CollectSvMetrics.size_vs_af_distrib, []])
  Array[File] joint_distrib_shards = flatten([short_joint_distribs, sv_joint_distribs])
  if ( length(joint_distrib_shards) > 0 ) {
    call QcTasks.SumCompressedDistribs as SumJointDistribs {
      input:
        distrib_tsvs = joint_distrib_shards,
        n_key_columns = 3,
        out_prefix = output_prefix + ".size_vs_af_distribution",
        g2c_analysis_docker = g2c_analysis_docker
    }
  }

  output {
    File? snv_sites_bed = ConcatSnvs.merged_file
    File? indel_sites_bed = ConcatIndels.merged_file
    File? sv_sites_bed = ConcatSvs.merged_file
    File? size_distrib = SumSizeDistribs.merged_distrib
    File? af_distrib = SumAfDistribs.merged_distrib
    File? size_vs_af_distrib = SumJointDistribs.merged_distrib
  }
}


# Slice & clean a gnomAD VCF to minimal info for site metric collection
task SliceVcf {
  input {
    File vcf
    File vcf_idx
    String interval

    String out_prefix

    Int disk_gb = 25
    String g2c_pipeline_docker
  }

  parameter_meta {
    vcf: {
      localization_optional: true
    }
  }

  command <<<
    set -eu -o pipefail

    # Parse intervals to set bounds on POS (important for avoiding double-counting SVs)
    n_int_parts=$( echo "~{interval}" | sed 's/:\|-/\n/g' | wc -l )
    if [ $n_int_parts -lt 3 ]; then
      min_pos=-1
      max_pos=1000000000
    else
      min_pos=$( echo "~{interval}" | sed 's/:\|-/\t/g' | cut -f2 )
      max_pos=$( echo "~{interval}" | sed 's/:\|-/\t/g' | cut -f3 )
    fi

    # Symlink vcf_idx to current working dir
    ln -s ~{vcf_idx} .

    # Ensure all necessary fields are defined in VCF header
    echo "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" > header.supp.vcf
    echo "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">" >> header.supp.vcf
    echo "##INFO=<ID=CN_NONREF_FREQ,Number=1,Type=Float,Description=\"CNV frequency\">" >> header.supp.vcf
    echo "##INFO=<ID=CN_NONREF_COUNT,Number=1,Type=Integer,Description=\"Nondip count.\">" >> header.supp.vcf
    echo "##INFO=<ID=AC_Het,Number=A,Type=Integer,Description=\"Heterozygous allele counts\">" >> header.supp.vcf
    echo "##INFO=<ID=AC_Hom,Number=A,Type=Integer,Description=\"Homozygous allele counts\">" >> header.supp.vcf
    echo "##INFO=<ID=AC_Hemi,Number=A,Type=Integer,Description=\"Hemizygous allele counts\">" >> header.supp.vcf
    echo "##INFO=<ID=HWE,Number=A,Type=Float,Description=\"HWE test\">" >> header.supp.vcf
    echo "##FILTER=<ID=MULTIALLELIC,Description=\"Multiallelic site\">" >> header.supp.vcf

    # Stream VCF to interval of interest and remove all INFO except for 
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
    bcftools view \
      --regions "~{interval}" \
      ~{vcf} \
    | bcftools annotate \
      -h header.supp.vcf \
    | bcftools view \
      --include '(FILTER = "PASS" & INFO/AC > 0) | INFO/SVTYPE = "CNV" | FILTER = "MULTIALLELIC"' \
    | awk -v min_pos="$min_pos" -v max_pos="$max_pos" \
      '{ if ($1 ~ "^#" || ($2 >= min_pos && $2 <= max_pos)) print }' \
    | bcftools annotate \
      -x "^INFO/END,INFO/SVTYPE,INFO/SVLEN,INFO/AN,INFO/AC,INFO/AF,INFO/CN_NONREF_COUNT,INFO/CN_NONREF_FREQ,INFO/AC_Het,INFO/AC_Hom,INFO/AC_Hemi,INFO/HWE" \
      -Oz -o "~{out_prefix}.cleaned.vcf.gz"
    tabix -p vcf "~{out_prefix}.cleaned.vcf.gz"

  >>>

  output {
    File cleaned_vcf = "~{out_prefix}.cleaned.vcf.gz"
    File cleaned_vcf_idx = "~{out_prefix}.cleaned.vcf.gz.tbi"
  }

  runtime {
    docker: g2c_pipeline_docker
    memory: "2 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    max_retries: 2
  }
}


task SetDefaultMinAf {
  input {
    Int? snv_n_samples
    Int? sv_n_samples
    String docker
  }

  command <<<
    python3 << 'EOF'
import math

# Always interpolate as *strings*, then decode in Python
snv_raw = "~{default="NA" snv_n_samples}"
sv_raw  = "~{default="NA" sv_n_samples}"
snv = None if snv_raw == "NA" else int(snv_raw)
sv  = None if sv_raw == "NA" else int(sv_raw)

vals = []
for x in [snv, sv]:
    if x is not None and x > 0:
        vals.append(1.0 / (2.0 * x))

if vals:
    result = max(vals)
else:
    result = 0.0  # or some other sentinel default

with open("out.txt", "w") as f:
    f.write(str(result) + "\n")
EOF
  >>>

  output {
    Float default_min_af = read_float("out.txt")
  }

  runtime {
    docker: docker
    memory: "2 GB"
    cpu: 1
    disks: "local-disk 20 HDD"
    preemptible: 3
    max_retries: 1
  }
}

