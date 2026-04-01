# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Preprocess 1kGP VCFs across technologies to match G2C QC expectations


version 1.0


import "https://raw.githubusercontent.com/collins-genomics/pancan_germline_wgs/refs/heads/posthoc_qc/wdl/Utilities.wdl" as Utils
import "https://raw.githubusercontent.com/collins-genomics/pancan_germline_wgs/refs/heads/posthoc_qc/wdl/vcf-qc/QcTasks.wdl" as QcTasks


workflow Preprocess1kgpVcfs {
  input {
    String srwgs_snv_vcf_prefix = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_"
    String srwgs_snv_vcf_suffix = ".recalibrated_variants.vcf.gz"

    String srwgs_sv_vcf_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz"
    String srwgs_sv_tbi_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/1KGP_3202.gatksv_svtools_novelins.freeze_V3.wAF.vcf.gz.tbi"

    String lrwgs_snv_vcf_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/vcf_with_unphased/variants_GRCh38_snv_snv_alt_with-unphased_HGSVC2024v1.0.vcf.gz"
    String lrwgs_snv_tbi_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/vcf_with_unphased/variants_GRCh38_snv_snv_alt_with-unphased_HGSVC2024v1.0.vcf.gz.tbi"

    String lrwgs_indel_vcf_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/vcf_with_unphased/variants_GRCh38_indel_insdel_alt_with-unphased_HGSVC2024v1.0.vcf.gz"
    String lrwgs_indel_tbi_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/vcf_with_unphased/variants_GRCh38_indel_insdel_alt_with-unphased_HGSVC2024v1.0.vcf.gz.tbi"

    String lrwgs_cnv_vcf_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/vcf_with_unphased/variants_GRCh38_sv_insdel_sym_with-unphased_HGSVC2024v1.0.vcf.gz"
    String lrwgs_cnv_tbi_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/vcf_with_unphased/variants_GRCh38_sv_insdel_sym_with-unphased_HGSVC2024v1.0.vcf.gz.tbi"

    String lrwgs_inv_vcf_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/vcf_with_unphased/variants_GRCh38_sv_inv_sym_with-unphased_HGSVC2024v1.0.vcf.gz"
    String lrwgs_inv_tbi_url = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/GRCh38/vcf_with_unphased/variants_GRCh38_sv_inv_sym_with-unphased_HGSVC2024v1.0.vcf.gz.tbi"

    File ref_fasta
    File ref_fasta_idx

    Array[File] srwgs_snv_scatter_intervals # Expects one file per chromosome in the same order as `contigs`
    Array[String] contigs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", 
                             "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                             "chr13", "chr14", "chr15", "chr16", "chr17", 
                             "chr18", "chr19", "chr20", "chr21", "chr22", 
                             "chrX",  "chrY"]
    
    String dest_bucket_uri_prefix = "gs://dfci-g2c-refs/hgsv/"

    String g2c_pipeline_docker
  }

  call QcTasks.MakeHeaderFiller {}

  scatter ( contig_info in zip(contigs, srwgs_snv_scatter_intervals) ) {

    String contig = contig_info.left
    File snv_scatter_interval = contig_info.right
    String vcf_ftp_url = srwgs_snv_vcf_prefix + contig + srwgs_snv_vcf_suffix
    String tbi_ftp_url = srwgs_snv_vcf_prefix + contig + srwgs_snv_vcf_suffix + ".tbi"

    call Utils.FtpDownload as DownloadSrwgsSnvVcf {
      input:
        ftp_url = vcf_ftp_url,
        output_name = basename(vcf_ftp_url),
        lftp_docker = g2c_pipeline_docker,
        disk_gb = 200
    }

    call Utils.FtpDownload as DownloadSrwgsSnvTbi {
      input:
        ftp_url = tbi_ftp_url,
        output_name = basename(tbi_ftp_url),
        lftp_docker = g2c_pipeline_docker,
        disk_gb = 200
    }

    File srwgs_snv_vcf = DownloadSrwgsSnvVcf.downloaded_file
    File srwgs_snv_tbi = DownloadSrwgsSnvTbi.downloaded_file

    call QcTasks.ParseIntervals as MakeSrwgsSnvIntervals {
      input:
        intervals_list = snv_scatter_interval,
        docker = "ubuntu:plucky-20251001"
    }
    
    scatter ( snv_interval_info in MakeSrwgsSnvIntervals.interval_info ) {

      String snv_vcf_shard_suffix = snv_interval_info[0]
      String snv_interval_coords = snv_interval_info[1]
      String snv_shard_fname = "srwgs.snv." + contig + "." + snv_vcf_shard_suffix + ".vcf.gz"

      call CurateSrwgsSnvs {
        input:
          vcf = srwgs_snv_vcf,
          vcf_tbi = srwgs_snv_tbi,
          out_vcf_fname = snv_shard_fname,
          interval = snv_interval_coords,
          ref_fasta = ref_fasta,
          ref_fasta_idx = ref_fasta_idx,
          supp_vcf_header = MakeHeaderFiller.supp_vcf_header,
          g2c_pipeline_docker = g2c_pipeline_docker
      }
    }

    call Utils.ConcatVcfs as ConcatSrwgsSnvs {
      input:
        vcfs = CurateSrwgsSnvs.vcf_out,
        vcf_idxs = CurateSrwgsSnvs.tbi_out,
        out_prefix = "1KGP.srWGS.snv_indel.cleaned." + contig,
        bcftools_concat_options = "--allow-overlaps --threads 2",
        bcftools_docker = g2c_pipeline_docker
    }

    call Utils.CopyGcpObjects as CopySrwgsSnvVcf {
      input:
        files_to_copy = [ConcatSrwgsSnvs.merged_vcf, ConcatSrwgsSnvs.merged_vcf_idx],
        destination = dest_bucket_uri_prefix + "dense_vcfs/srwgs/snv_indel/"
    }

    call CurateSrwgsSvs {
      input:
        vcf_url = srwgs_sv_vcf_url,
        tbi_url = srwgs_sv_tbi_url,
        dest_bucket = dest_bucket_uri_prefix + "dense_vcfs/srwgs/sv/",
        contig = contig,
        supp_vcf_header = MakeHeaderFiller.supp_vcf_header,
        g2c_pipeline_docker = g2c_pipeline_docker
    }

    call CurateLrwgsSnvs {
      input:
        snv_vcf_url = lrwgs_snv_vcf_url,
        snv_tbi_url = lrwgs_snv_tbi_url,
        indel_vcf_url = lrwgs_indel_vcf_url,
        indel_tbi_url = lrwgs_indel_tbi_url,
        dest_bucket = dest_bucket_uri_prefix + "dense_vcfs/lrwgs/snv_indel/",
        contig = contig,
        ref_fasta = ref_fasta,
        ref_fasta_idx = ref_fasta_idx,
        supp_vcf_header = MakeHeaderFiller.supp_vcf_header,
        g2c_pipeline_docker = g2c_pipeline_docker
    }

    call CurateLrwgsSvs {
      input:
        cnv_vcf_url = lrwgs_cnv_vcf_url,
        cnv_tbi_url = lrwgs_cnv_tbi_url,
        inv_vcf_url = lrwgs_inv_vcf_url,
        inv_tbi_url = lrwgs_inv_tbi_url,
        dest_bucket = dest_bucket_uri_prefix + "dense_vcfs/lrwgs/sv/",
        contig = contig,
        supp_vcf_header = MakeHeaderFiller.supp_vcf_header,
        g2c_pipeline_docker = g2c_pipeline_docker
    }
  }
}


task CurateSrwgsSnvs {
  input {
    File vcf
    File vcf_tbi
    String out_vcf_fname
    String interval

    File ref_fasta
    File ref_fasta_idx

    File supp_vcf_header

    String g2c_pipeline_docker
  }

  Int disk_gb = ceil(2 * size(vcf, "GB")) + 10
  String out_tbi_fname = out_vcf_fname + ".tbi"

  command <<<
    set -eu -o pipefail

    # Need to copy VCF index to same directory as VCF on VM
    # This is due to a weird interaction between how Cromwell localizes input files
    # and where bcftools searches for index files
    cp ~{vcf_tbi} ~{vcf}.tbi

    # Filter & clean VCF
    bcftools view --regions "~{interval}" ~{vcf} \
    | bcftools norm -m -any -f ~{ref_fasta} -c s --threads 2 \
    | bcftools view -f PASS -c 1 \
    | bcftools +fill-tags \
    | bcftools annotate -h ~{supp_vcf_header} \
    | bcftools annotate --threads 2 \
      --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
      -x "^INFO/END,INFO/SVTYPE,INFO/SVLEN,INFO/AN,INFO/AC,INFO/AF,INFO/CN_NONREF_COUNT,INFO/CN_NONREF_FREQ,INFO/AC_Het,INFO/AC_Hom,INFO/AC_Hemi,INFO/HWE,^FORMAT/GT,^FILTER/PASS" \
      -Oz -o "~{out_vcf_fname}"
    tabix -f -p vcf "~{out_vcf_fname}"
  >>>

  output {
    File vcf_out = "~{out_vcf_fname}"
    File tbi_out = "~{out_tbi_fname}"
  }

  runtime {
    docker: g2c_pipeline_docker
    memory: "4 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    max_retries: 1
  }
}


task CurateSrwgsSvs {
  input {
    String vcf_url
    String tbi_url
    String dest_bucket
    String contig

    File supp_vcf_header

    String g2c_pipeline_docker
    Int disk_gb = 30
  }

  String in_vcf_fname = basename(vcf_url)

  String out_vcf_fname = "1KGP.srWGS.sv.cleaned." + contig + ".vcf.gz"
  String out_tbi_fname = out_vcf_fname + ".tbi"

  command <<<
    set -eu -o pipefail

    # Download VCF and index
    wget ~{vcf_url}
    wget ~{tbi_url}

    # Filter & clean VCF
    bcftools view --regions "~{contig}" -f PASS,MULTIALLELIC --include 'INFO/SOURCE="gatksv"' "~{in_vcf_fname}" \
    | bcftools +fill-tags \
    | bcftools annotate -h ~{supp_vcf_header} \
    | bcftools annotate --threads 2 \
      -x "^INFO/END,INFO/SVTYPE,INFO/SVLEN,INFO/AN,INFO/AC,INFO/AF,INFO/CN_NONREF_COUNT,INFO/CN_NONREF_FREQ,INFO/AC_Het,INFO/AC_Hom,INFO/AC_Hemi,INFO/HWE,^FORMAT/GT,FORMAT/RD_CN,^FILTER/PASS,FILTER/MULTIALLELIC" \
    | /opt/pancan_germline_wgs/scripts/data_management/external_data_curation/postprocess_hgsvc_vcfs.py --minimal-sv --no-gt stdin "~{out_vcf_fname}"
    tabix -f -p vcf "~{out_vcf_fname}"

    # Stage cleaned VCF to desired output bucket
    gsutil -m cp \
      "~{out_vcf_fname}" \
      "~{out_tbi_fname}" \
      ~{dest_bucket}
  >>>

  output {}

  runtime {
    docker: g2c_pipeline_docker
    memory: "4 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    max_retries: 1
  }
}


task CurateLrwgsSnvs {
  input {
    String snv_vcf_url
    String snv_tbi_url
    String indel_vcf_url
    String indel_tbi_url
    String dest_bucket
    String contig

    File ref_fasta
    File ref_fasta_idx

    File supp_vcf_header

    String g2c_pipeline_docker
    Int disk_gb = 20
  }

  String in_snv_fname = basename(snv_vcf_url)
  String in_indel_fname = basename(indel_vcf_url)

  String out_vcf_fname = "1KGP.lrWGS.snv_indel.cleaned." + contig + ".vcf.gz"
  String out_tbi_fname = out_vcf_fname + ".tbi"

  command <<<
    set -eu -o pipefail

    # Download VCFs and indexes
    wget ~{snv_vcf_url}
    wget ~{snv_tbi_url}
    wget ~{indel_vcf_url}
    wget ~{indel_tbi_url}

    # Merge, filter, and clean VCFs
    bcftools concat --allow-overlaps --regions ~{contig} ~{in_snv_fname} ~{in_indel_fname} \
    | bcftools norm -m -any -f ~{ref_fasta} -c s --threads 2 \
    | bcftools annotate -h ~{supp_vcf_header} -x "^INFO/VARTYPE,^FORMAT/GT" \
    | /opt/pancan_germline_wgs/scripts/data_management/external_data_curation/postprocess_hgsvc_vcfs.py stdin stdout \
    | bcftools view -c 1 \
    | bcftools +fill-tags \
    | bcftools annotate --threads 2 \
      --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
      -x "^INFO/END,INFO/SVTYPE,INFO/SVLEN,INFO/AN,INFO/AC,INFO/AF,INFO/CN_NONREF_COUNT,INFO/CN_NONREF_FREQ,INFO/AC_Het,INFO/AC_Hom,INFO/AC_Hemi,INFO/HWE,^FORMAT/GT" \
      -Oz -o "~{out_vcf_fname}"
    tabix -f -p vcf "~{out_vcf_fname}"

    # Stage cleaned VCF to desired output bucket
    gsutil -m cp \
      "~{out_vcf_fname}" \
      "~{out_tbi_fname}" \
      ~{dest_bucket}
  >>>

  output {}

  runtime {
    docker: g2c_pipeline_docker
    memory: "4 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    max_retries: 1
  }
}


task CurateLrwgsSvs {
  input {
    String cnv_vcf_url
    String cnv_tbi_url
    String inv_vcf_url
    String inv_tbi_url
    String dest_bucket
    String contig

    File supp_vcf_header

    String g2c_pipeline_docker
    Int disk_gb = 30
  }

  String in_cnv_fname = basename(cnv_vcf_url)
  String in_inv_fname = basename(inv_vcf_url)

  String out_vcf_fname = "1KGP.lrWGS.sv.cleaned." + contig + ".vcf.gz"
  String out_tbi_fname = out_vcf_fname + ".tbi"

  command <<<
    set -eu -o pipefail

    # Download VCFs and indexes
    wget ~{cnv_vcf_url}
    wget ~{cnv_tbi_url}
    wget ~{inv_vcf_url}
    wget ~{inv_tbi_url}

    # Merge, filter, and clean VCFs
    bcftools concat --allow-overlaps --regions ~{contig} ~{in_cnv_fname} ~{in_inv_fname} \
    | bcftools annotate -h ~{supp_vcf_header} -x "^INFO/SVTYPE,INFO/SVLEN,^FORMAT/GT" \
    | /opt/pancan_germline_wgs/scripts/data_management/external_data_curation/postprocess_hgsvc_vcfs.py stdin stdout --sv \
    | bcftools +fill-tags \
    | bcftools view -c 1 \
    | bcftools annotate --threads 2 \
      -x "^INFO/END,INFO/SVTYPE,INFO/SVLEN,INFO/AN,INFO/AC,INFO/AF,INFO/CN_NONREF_COUNT,INFO/CN_NONREF_FREQ,INFO/AC_Het,INFO/AC_Hom,INFO/AC_Hemi,INFO/HWE,^FORMAT/GT" \
      -Oz -o "~{out_vcf_fname}"
    tabix -f -p vcf "~{out_vcf_fname}"

    # Stage cleaned VCF to desired output bucket
    gsutil -m cp \
      "~{out_vcf_fname}" \
      "~{out_tbi_fname}" \
      ~{dest_bucket}
  >>>

  output {}

  runtime {
    docker: g2c_pipeline_docker
    memory: "4 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    max_retries: 1
  }
}
