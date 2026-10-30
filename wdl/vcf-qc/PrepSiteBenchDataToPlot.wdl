# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Helper subroutine to preprocess site benchmarking data prior to QC plotting
# Not really intended to be called as a stand-alone workflow, 
# but is called from within PlotVcfQcMetrics.wdl


version 1.0


import "QcTasks.wdl" as QcTasks


workflow PrepSiteBenchDataToPlot {
  input {
    String dataset_name
    String dataset_prefix
    Array[String] interval_set_names = ["all"]

    Array[Array[File?]] common_snv_ppv_beds = [[]]
    Array[Array[File?]] common_indel_ppv_beds = [[]]
    Array[Array[File?]] common_sv_ppv_beds = [[]]

    Array[Array[File?]] ppv_by_freqs = [[]]
    Array[Array[File?]] sensitivity_by_freqs = [[]]
    Int n_key_columns = 3

    String bcftools_docker
    String g2c_analysis_docker
  }

  Int n_interval_sets = length(interval_set_names)
  Int n_snv_beds_any = length(select_all(flatten(common_snv_ppv_beds)))
  Int n_indel_beds_any = length(select_all(flatten(common_indel_ppv_beds)))
  Int n_sv_beds_any = length(select_all(flatten(common_sv_ppv_beds)))

  # Process each set of evaluation intervals in parallel
  scatter ( sbi in range(n_interval_sets) ) {

    String bi_name = interval_set_names[sbi]
    String sbi_prefix = dataset_prefix + "." + bi_name

    Int n_snv_beds = if length(common_snv_ppv_beds) > sbi then length(select_all(common_snv_ppv_beds[sbi])) else 0
    Int n_indel_beds = if length(common_indel_ppv_beds) > sbi then length(select_all(common_indel_ppv_beds[sbi])) else 0
    Int n_sv_beds = if length(common_sv_ppv_beds) > sbi then length(select_all(common_sv_ppv_beds[sbi])) else 0
    Int n_ppv_tsvs = if length(ppv_by_freqs) > sbi then length(select_all(ppv_by_freqs[sbi])) else 0
    Int n_sens_tsvs = if length(sensitivity_by_freqs) > sbi then length(select_all(sensitivity_by_freqs[sbi])) else 0
      
    # Collapse SNV BEDs
    if (n_snv_beds > 0) {
      Array[File] sb_common_snv_preflat = select_all(common_snv_ppv_beds[sbi])
      if ( n_snv_beds > 1 ) {
        call QcTasks.ConcatTextFiles as CollapseSiteBenchSnvs {
          input:
            shards = sb_common_snv_preflat,
            concat_command = "zcat",
            sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
            compression_command = "bgzip -c",
            input_has_header = true,
            output_filename = sbi_prefix + ".common_snvs.bed.gz",
            docker = bcftools_docker
        }
      }
      File sb_common_snv_flat = select_first([CollapseSiteBenchSnvs.merged_file,
                                              sb_common_snv_preflat[0]])
    }
      
    # Collapse indel BEDs
    if (n_indel_beds > 0) {
      Array[File] sb_common_indel_preflat = select_all(common_indel_ppv_beds[sbi])
      if ( n_indel_beds > 1 ) {
        call QcTasks.ConcatTextFiles as CollapseSiteBenchIndels {
          input:
            shards = sb_common_indel_preflat,
            concat_command = "zcat",
            sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
            compression_command = "bgzip -c",
            input_has_header = true,
            output_filename = sbi_prefix + ".common_indels.bed.gz",
            docker = bcftools_docker
        }
      }
      File sb_common_indel_flat = select_first([CollapseSiteBenchIndels.merged_file,
                                                sb_common_indel_preflat[0]])
    }
      
    # Collapse SV BEDs
    if (n_sv_beds > 0) {
      Array[File] sb_common_sv_preflat = select_all(common_sv_ppv_beds[sbi])
      if ( n_sv_beds > 1 ) {
        call QcTasks.ConcatTextFiles as CollapseSiteBenchSvs {
          input:
            shards = sb_common_sv_preflat,
            concat_command = "zcat",
            sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
            compression_command = "bgzip -c",
            input_has_header = true,
            output_filename = sbi_prefix + ".common_svs.bed.gz",
            docker = bcftools_docker
        }
      }
      File sb_common_sv_flat = select_first([CollapseSiteBenchSvs.merged_file,
                                             sb_common_sv_preflat[0]])
    }

    # Collapse compressed PPV tsvs
    if (n_ppv_tsvs > 0) {
      call QcTasks.SumCompressedDistribs as SumSiteBenchPpvByAf {
        input:
          distrib_tsvs = select_all(ppv_by_freqs[sbi]),
          n_key_columns = n_key_columns,
          out_prefix = sbi_prefix + ".ppv_by_freq",
          g2c_analysis_docker = g2c_analysis_docker
      }
    }

    # # Collapse compressed sensitivity tsvs
    if (n_sens_tsvs > 0) {
      call QcTasks.SumCompressedDistribs as SumSiteBenchSensByAf {
        input:
          distrib_tsvs = select_all(sensitivity_by_freqs[sbi]),
          n_key_columns = n_key_columns,
          out_prefix = sbi_prefix + ".sensitivity_by_freq",
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
  }

  # Further collapse SNVs across interval sets
  if ( n_snv_beds_any > 0 ) {
    Array[File] sb_common_snv_flat_array = select_all(sb_common_snv_flat)
    if ( length(sb_common_snv_flat_array) > 1 ) {
      call QcTasks.ConcatTextFiles as CollapseSiteBenchSnvsUnion {
        input:
          shards = sb_common_snv_flat_array,
          concat_command = "zcat",
          sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
          compression_command = "bgzip -c",
          input_has_header = true,
          output_filename = dataset_prefix + ".merged.common_snvs.bed.gz",
          docker = bcftools_docker
      }
    }
    File sb_common_snv_flat_union = select_first([CollapseSiteBenchSnvsUnion.merged_file,
                                                  sb_common_snv_flat_array[0]])
  }

  # Further collapse indels across interval sets
  if ( n_indel_beds_any > 0 ) {
    Array[File] sb_common_indel_flat_array = select_all(sb_common_indel_flat)
    if ( length(sb_common_indel_flat_array) > 1 ) {
      call QcTasks.ConcatTextFiles as CollapseSiteBenchIndelsUnion {
        input:
          shards = sb_common_indel_flat_array,
          concat_command = "zcat",
          sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
          compression_command = "bgzip -c",
          input_has_header = true,
          output_filename = dataset_prefix + ".merged.common_indels.bed.gz",
          docker = bcftools_docker
      }
    }
    File sb_common_indel_flat_union = select_first([CollapseSiteBenchIndelsUnion.merged_file,
                                                  sb_common_indel_flat_array[0]])
  }

  # Further collapse SVs across interval sets
  if ( n_sv_beds_any > 0 ) {
    Array[File] sb_common_sv_flat_array = select_all(sb_common_sv_flat)
    if ( length(sb_common_sv_flat_array) > 1 ) {
      call QcTasks.ConcatTextFiles as CollapseSiteBenchSvsUnion {
        input:
          shards = sb_common_sv_flat_array,
          concat_command = "zcat",
          sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
          compression_command = "bgzip -c",
          input_has_header = true,
          output_filename = dataset_prefix + ".merged.common_svs.bed.gz",
          docker = bcftools_docker
      }
    }
    File sb_common_sv_flat_union = select_first([CollapseSiteBenchSvsUnion.merged_file,
                                                  sb_common_sv_flat_array[0]])
  }

  # For convenience, dump all output arrays into a .json that can be parsed 
  # by downstream plotting tasks
  call MakeOutputMenu {
    input:
      snv_ppv_bed_uris = sb_common_snv_flat,
      indel_ppv_bed_uris = sb_common_indel_flat,
      sv_ppv_bed_uris = sb_common_sv_flat,
      snv_ppv_union_uri = sb_common_snv_flat_union,
      indel_ppv_union_uri = sb_common_indel_flat_union,
      sv_ppv_union_uri = sb_common_sv_flat_union,
      ppv_by_af_uris = SumSiteBenchPpvByAf.merged_distrib,
      sens_by_af_uris = SumSiteBenchSensByAf.merged_distrib,
      out_json = dataset_prefix + ".site_benchmarking_plot_inputs.json",
      python_docker = g2c_analysis_docker
  }

  output {
    File plot_files_json = MakeOutputMenu.json_out

    Array[File?] snv_ppv_beds = sb_common_snv_flat
    Array[File?] indel_ppv_beds = sb_common_indel_flat
    Array[File?] sv_ppv_beds = sb_common_sv_flat

    File? snv_ppv_union = sb_common_snv_flat_union
    File? indel_ppv_union = sb_common_indel_flat_union
    File? sv_ppv_union_union = sb_common_sv_flat_union
  
    Array[File?] ppv_by_af = SumSiteBenchPpvByAf.merged_distrib
    Array[File?] sens_by_af = SumSiteBenchSensByAf.merged_distrib
  }
}


# Dump the URIs of all prepared site benchmarking files to a .json
task MakeOutputMenu {
  input {
    Array[String?] snv_ppv_bed_uris
    Array[String?] indel_ppv_bed_uris
    Array[String?] sv_ppv_bed_uris

    String? snv_ppv_union_uri
    String? indel_ppv_union_uri
    String? sv_ppv_union_uri

    Array[String?] ppv_by_af_uris
    Array[String?] sens_by_af_uris

    String out_json
    String python_docker = "python:3.9-slim"
  }

  File snv_ppv_json = write_json(select_first([snv_ppv_bed_uris, []]))
  File indel_ppv_json = write_json(select_first([indel_ppv_bed_uris, []]))
  File sv_ppv_json = write_json(select_first([sv_ppv_bed_uris, []]))
  File ppv_by_af_json = write_json(select_first([ppv_by_af_uris, []]))
  File sens_by_af_json = write_json(select_first([sens_by_af_uris, []]))
  String snv_ppv_union_uri_use = if defined(snv_ppv_union_uri)
                                 then "\"" + select_first([snv_ppv_union_uri, ""]) + "\""
                                 else "null"
  String indel_ppv_union_uri_use = if defined(indel_ppv_union_uri)
                                   then "\"" + select_first([indel_ppv_union_uri, ""]) + "\""
                                   else "null"
  String sv_ppv_union_uri_use = if defined(sv_ppv_union_uri)
                                then "\"" + select_first([sv_ppv_union_uri, ""]) + "\""
                                else "null"

  command <<<
    set -eu -o pipefail

    cat > ~{out_json} <<EOF
{
  "snv_ppv_beds": $( cat ~{snv_ppv_json} ),
  "indel_ppv_beds": $( cat ~{indel_ppv_json} ),
  "sv_ppv_beds": $( cat ~{sv_ppv_json} ),
  "snv_ppv_union": ~{snv_ppv_union_uri_use},
  "indel_ppv_union": ~{indel_ppv_union_uri_use},
  "sv_ppv_union": ~{sv_ppv_union_uri_use},
  "ppv_by_af": $( cat ~{ppv_by_af_json} ),
  "sens_by_af": $( cat ~{sens_by_af_json} )
}
EOF

    # For ease of debugging, print the output .json to stdout
    cat ~{out_json} | python -m json.tool
  >>>

  output {
    File json_out = "~{out_json}"
  }

  runtime {
    docker: python_docker
    memory: "1.7 GB"
    cpu: 1
    disks: "local-disk 20 HDD"
    preemptible: 3
    max_retries: 1
  }
}
