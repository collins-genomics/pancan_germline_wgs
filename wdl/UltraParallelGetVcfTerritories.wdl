# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Bulk chromosomal territory computation for a very large VCF array

# This is a helper wrapper around GetGenomeTerritoryPerVcf.wdl
# designed to handle very large VCF shard arrays.
# Please refer to that WDL for details on the method itself.


version 1.0


import "Utilities.wdl" as Utils
import "GetGenomeTerritoryPerVcf.wdl" as GetTerritories


workflow UltraParallelGetVcfTerritories {
	input {
		File vcf_uri_list         # List of GCS URIs for all VCFs to be processed. Tabix indexes for each must exist in the same bucket
    File genome_file          # BEDTools-style genome file
    Int vcfs_per_chunk = 100  # We have experienced Cromwell input choking at around ~1,000. Recommended to keep this parameter <500

    String output_prefix

    String g2c_analysis_docker
	}

  # Shard vcf_uri_list
  call Utils.ShardTextFile as ShardVcfList {
    input:
      input_file = vcf_uri_list,
      lines_per_split = vcfs_per_chunk,
      out_prefix = output_prefix + ".",
      g2c_analysis_docker = g2c_analysis_docker
  }
  Int n_chunks = length(ShardVcfList.shards)

  # Scatter over sharded VCF lists
  scatter ( i in range(n_chunks) ) {
    # Extract URIs from sharded file as array of strings and indexes
    call Utils.ExtractVcfArrays {
      input:
        vcf_list = ShardVcfList.shards[i],
        linux_docker = g2c_analysis_docker
    }

    # Get territories for all VCFs in this chunk
    call GetTerritories.GetGenomeTerritoryPerVcf {
      input:
        vcfs = ExtractVcfArrays.vcf_uris,
        vcf_idxs = ExtractVcfArrays.vcf_tbi_uris,
        compute_density = false,
        genome_file = genome_file,
        output_prefix = "~{output_prefix}.chunk_~{i}",
        g2c_analysis_docker = g2c_analysis_docker
    }
  }

  # Further combine territory density maps across chunks
  call Utils.CalcBedDensity as CalcDensity {
    input:
      beds = GetGenomeTerritoryPerVcf.territories,
      genome_file = select_first(select_all([genome_file])),
      output_prefix = output_prefix,
      bedtools_docker = g2c_analysis_docker
  }

  output {
    File territory_density = CalcDensity.density_bed
  }
}
