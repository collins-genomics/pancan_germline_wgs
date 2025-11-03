# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Generic WDL tasks used for VCF quality control

# Note that due to discrepancies in how cromshell and Dockstore handle
# resolve relative imports, some functions here were duplicated 
# from Utilities.wdl in the parent wdl/ subdirectory; this was 
# deemed the "least bad" solution to this problem given our constraints


version 1.0


# Compare genotypes between callsets for a list of samples
task BenchmarkGenotypes {
  input {
    File variant_id_map
    File sample_id_map
    Boolean invert_sample_map = false
    File source_site_metrics
    String output_prefix

    File source_gt_tarball
    File target_gt_tarball

    Boolean report_by_gt = false
    Float common_af_cutoff = 0.01

    Float mem_gb = 7.5
    Int n_cpu = 4
    
    String g2c_analysis_docker
  }

  parameter_meta {
    source_gt_tarball: {
      localization_optional: true
    }
    target_gt_tarball: {
      localization_optional: true
    }
  }

  String gt_report_cmd = if report_by_gt then "--report-by-genotype" else ""
  String invert_sid_cmd = if invert_sample_map then "--invert-sid" else ""
  Int disk_gb = ceil(5 * size([source_gt_tarball, target_gt_tarball, source_site_metrics], "GB")) + 20

  command <<<
    set -eu -o pipefail

    # Prep variant ID lists & metric files
    touch source.vids.list target.vids.list
    zcat ~{variant_id_map} | cut -f1 | sort | uniq > source.vids.list || true
    zcat ~{variant_id_map} | cut -f2 | fgrep -xv "NA" | sort | uniq > target.vids.list || true

    # If there are no source and target variants, there's no need to run the rest of this task
    if [ $( cat source.vids.list | wc -l ) -eq 0 ] && \
       [ $( cat target.vids.list | wc -l ) -eq 0 ]; then

      if ~{report_by_gt}; then
        echo -e "#sample\tclass\tsubclass\tfreq_bin\tgenotype\tno_match\tcarrier_match\tgt_match" \
        > ~{output_prefix}.gt_comparison.distrib.tsv
      else
        echo -e "#sample\tclass\tsubclass\tfreq_bin\tzygosity\tno_match\tcarrier_match\tgt_match" \
        > ~{output_prefix}.gt_comparison.distrib.tsv
      fi

    else

      # Prep variant metric files
      zcat ~{source_site_metrics} | sed -n '1p' | fgrep "#" > site_metrics.header
      zcat ~{source_site_metrics} | fgrep -wf source.vids.list \
      | cat site_metrics.header - | bgzip -c \
      > source.metrics.bed.gz
      rm ~{source_site_metrics}

      # Prep sample lists
      if ~{invert_sample_map}; then
        awk -v FS="\t" -v OFS="\t" '{ print $2, $1 }' ~{sample_id_map} > sample.map.tsv
      else
        cp ~{sample_id_map} sample.map.tsv
      fi
      cut -f1 sample.map.tsv | sort > source.samples.list
      cut -f2 sample.map.tsv | sort > target.samples.list

      # Unpack source GTs and subset to samples of interest
      mkdir source_gts_raw
      gsutil -m cp ~{source_gt_tarball} source_gt_tarball.tar.gz
      tar -xzvf source_gt_tarball.tar.gz -C source_gts_raw/
      mkdir source_gts/
      while read sid; do
        find source_gts_raw/ -name "$sid.gt.tsv.gz" \
        | xargs -I {} zcat {} | sort -k1 \
        | join -t$'\t' - source.vids.list \
        | gzip -c > source_gts/$sid.gt.sub.tsv.gz
      done < source.samples.list
      rm -rf source_gts_raw source_gt_tarball.tar.gz
      echo "Contents of source_gts:"
      ls -lh source_gts/

      # Unpack target GTs and subset to samples of interest
      mkdir target_gts_raw
      gsutil -m cp ~{target_gt_tarball} target_gt_tarball.tar.gz
      tar -xzvf target_gt_tarball.tar.gz -C target_gts_raw/
      mkdir target_gts/
      while read sid; do
        find target_gts_raw/ -name "$sid.gt.tsv.gz" \
        | xargs -I {} zcat {} | sort -k1 \
        | join -t$'\t' - target.vids.list \
        | gzip -c > target_gts/$sid.gt.sub.tsv.gz
      done < target.samples.list
      rm -rf target_gts_raw target_gt_tarball.tar.gz
      echo "Contents of target_gts:"
      ls -lh target_gts/

      # Benchmark genotypes
      echo "Now benchmarking..."
      /opt/pancan_germline_wgs/scripts/qc/vcf_qc/compare_genotypes.R \
        --variant-map ~{variant_id_map} \
        --source-site-metrics source.metrics.bed.gz \
        --sample-map sample.map.tsv \
        --source-gt-dir source_gts/ \
        --target-gt-dir target_gts/ \
        --gt-tsv-suffix ".gt.sub.tsv.gz" \
        ~{gt_report_cmd} \
        ~{invert_sid_cmd} \
        --common-af ~{common_af_cutoff} \
        --out-prefix ~{output_prefix}

    fi

    gzip -f ~{output_prefix}.gt_comparison.distrib.tsv
  >>>

  output {
    File gt_bench_distrib = "~{output_prefix}.gt_comparison.distrib.tsv.gz"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "~{mem_gb} GB"
    cpu: n_cpu
    maxRetries: 1
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}


task CollectSampleGenotypeMetrics {
  input {
    File vcf
    File vcf_idx

    File? site_metrics
    Float? common_af_cutoff

    String g2c_analysis_docker
  }

  String out_base = basename(vcf, ".vcf.gz")
  String gt_outfile = out_base + ".genotypes.tsv.gz"
  String distrib_outfile = out_base + ".gt_distrib.tsv.gz"

  String site_metrics_basename = if defined(site_metrics) then basename(select_first([site_metrics])) else ""
  String distrib_cmd = if defined(site_metrics) then "--site-metrics ~{site_metrics_basename} --distrib-out ~{out_base}.gt_distrib.tsv.gz" else ""
  String common_cmd = if defined(site_metrics) && defined(common_af_cutoff) then "--common-af ~{common_af_cutoff}" else ""

  Int disk_gb = ceil(1.2 * size(vcf, "GB")) + 5

  command <<<
    set -eu -o pipefail

    if ~{defined(site_metrics)}; then
      mv ~{select_first([site_metrics])} ./
    fi

    bcftools query -i 'GT="alt"' -f '[%CHROM\t%POS\t%REF\t%ALT\t%INFO/SVLEN\t%SAMPLE\t%GT\t%RD_CN\n]' ~{vcf} \
    | /opt/pancan_germline_wgs/scripts/qc/vcf_qc/clean_sample_genotypes.py ~{distrib_cmd} ~{common_cmd} \
    | gzip -c \
    > ~{gt_outfile}
  >>>

  output {
    File genotypes_tsv = gt_outfile
    File? compressed_gt_distrib = distrib_outfile
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}


task CollectSiteMetrics {
  input {
    File vcf
    File vcf_idx

    Int n_samples

    Float? min_af_bin
    Float? common_af_cutoff

    String g2c_analysis_docker
  }

  String out_prefix = basename(vcf, ".vcf.gz")
  Int disk_gb = ceil(2 * size(vcf, "GB")) + 5

  String min_af_cmd = if defined(min_af_bin) then "--min-af-bin ~{min_af_bin}" else ""
  String common_cmd = if defined(common_af_cutoff) then "--common-af ~{common_af_cutoff}" else ""

  command <<<
    set -eu -o pipefail

    # Ensure all necessary fields are defined in VCF header
    echo "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">" > header.supp.vcf
    echo "##INFO=<ID=CN_NONREF_FREQ,Number=1,Type=Float,Description=\"CNV frequency\">" >> header.supp.vcf
    echo "##INFO=<ID=CN_NONREF_COUNT,Number=1,Type=Integer,Description=\"Nondip count.\">" >> header.supp.vcf
    echo "##INFO=<ID=HWE,Number=A,Type=Float,Description=\"HWE test\">" >> header.supp.vcf
    echo "##INFO=<ID=ExcHet,Number=A,Type=Float,Description=\"ExcHet test\">" >> header.supp.vcf

    # Collect stats and split into SNV, indel, and SV files
    bcftools annotate -h header.supp.vcf ~{vcf} \
    | bcftools query \
      -f '%CHROM\t%POS\t%END\t%REF\t%ALT\t%INFO/SVLEN\t%INFO/AN\t%INFO/AC\t%INFO/AF\t%INFO/CN_NONREF_COUNT\t%INFO/CN_NONREF_FREQ\t%INFO/AC_Het\t%INFO/AC_Hom\t%INFO/AC_Hemi\t%INFO/HWE\n' \
    | /opt/pancan_germline_wgs/scripts/qc/vcf_qc/clean_site_metrics.py \
      ~{min_af_cmd} \
      ~{common_cmd} \
      -o ~{out_prefix} \
      --gzip \
      -N ~{n_samples}

    # Concatenate all site metrics for downstream 
    # compatability with CollectSampleGenotypeMetrics
    find ./ -name "~{out_prefix}.*.sites.bed.gz" > site_beds.list || true
    zcat $( head -n1 site_beds.list ) | head -n1 > site_metrics.header || true
    while read site_file; do
      zcat $site_file | fgrep -v "#" || true
    done < site_beds.list \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V -k6,6V \
    | uniq \
    | cat site_metrics.header - \
    | bgzip -c \
    > ~{out_prefix}.all.sites.bed.gz || true
    tabix -p bed -f ~{out_prefix}.all.sites.bed.gz

    # Concatenate all common variant IDs for downstrea
    # compatability with LD-based analyses
    find ./ -name "~{out_prefix}.*.sites.common.bed.gz" \
    | xargs -I {} zcat {} \
    | grep -ve '^#' | cut -f4 \
    | sort -V | uniq \
    > "~{out_prefix}.common_vids.list" || true
  >>>

  output {
    File all_sites = out_prefix + ".all.sites.bed.gz"
    File all_sites_idx = out_prefix + ".all.sites.bed.gz.tbi"
    File? snv_sites = out_prefix + ".snv.sites.bed.gz"
    File? snv_sites_idx = out_prefix + ".snv.sites.bed.gz.tbi"
    File? indel_sites = out_prefix + ".indel.sites.bed.gz"
    File? indel_sites_idx = out_prefix + ".indel.sites.bed.gz.tbi"
    File? sv_sites = out_prefix + ".sv.sites.bed.gz"
    File? sv_sites_idx = out_prefix + ".sv.sites.bed.gz.tbi"
    File? common_snv_sites = out_prefix + ".snv.sites.common.bed.gz"
    File? common_snv_sites_idx = out_prefix + ".snv.sites.common.bed.gz.tbi"
    File? common_indel_sites = out_prefix + ".indel.sites.common.bed.gz"
    File? common_indel_sites_idx = out_prefix + ".indel.sites.common.bed.gz.tbi"
    File? common_sv_sites = out_prefix + ".sv.sites.common.bed.gz"
    File? common_sv_sites_idx = out_prefix + ".sv.sites.common.bed.gz.tbi"
    File common_vids_list = out_prefix + ".common_vids.list"
    File size_distrib = out_prefix + ".size_distrib.tsv.gz"
    File af_distrib = out_prefix + ".af_distrib.tsv.gz"
    File size_vs_af_distrib = out_prefix + ".size_vs_af_distrib.tsv.gz"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "2 GB"
    cpu: 1
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
    maxRetries: 3
  }
}


# Concatenate an array of sample genotype .tsvs and split into one file per sample
task ConcatGenotypeTsvs {
  input {
    Array[File] tsvs
    String output_prefix
    String g2c_analysis_docker
  }

  Int disk_gb = ceil(2 * size(tsvs, "GB")) + 10
  String outdir = output_prefix + "_sample_genotypes"

  command <<<
    set -eu -o pipefail

    mkdir ~{outdir}

    # Split each input tsv by appending to sample-specific files in outdir/
    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/split_and_append_gts.py \
      ~{write_lines(tsvs)} \
      ~{outdir}

    # Compress outdir as tarball for easier GCP/Cromwell IO
    tar -czvf ~{outdir}.tar.gz ~{outdir}
  >>>

  output {
    File genotypes_tarball = "~{outdir}.tar.gz"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}


task ConcatTextFiles {
  input {
    Array[File] shards
    String concat_command = "cat"
    String? sort_command
    String? compression_command
    Boolean input_has_header = false
    String output_filename

    Float mem_gb = 1.75
    Int n_cpu = 1
    Int? disk_gb
    String docker
  }

  Int disk_gb_use = select_first([disk_gb, ceil(2 * size(shards, "GB")) + 10])
  String sort = if defined(sort_command) then " | " + select_first([sort_command, ""]) else ""
  String compress = if defined(compression_command) then " | " + select_first([compression_command, ""]) else ""
  String posthoc_cmds = if input_has_header then sort + " | fgrep -xvf header.txt | cat header.txt - " + compress else sort + compress

  command <<<
    set -euo pipefail

    if [ "~{input_has_header}" == "true" ]; then
      ~{concat_command} ~{shards[0]} \
      | head -n1 > header.txt || true
    else
      touch header.txt
    fi

    cat ~{write_lines(shards)} | xargs -I {} ~{concat_command} {} ~{posthoc_cmds} > ~{output_filename} || true
  >>>

  output {
    File merged_file = "~{output_filename}"
  }

  runtime {
    docker: docker
    memory: "~{mem_gb} GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb_use + " HDD"
    preemptible: 3
  }
}


# Duplicated from Utilities.wdl
task ConcatVcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    String out_prefix

    String bcftools_concat_options = ""

    Float mem_gb = 3.5
    Int cpu_cores = 2
    Int? disk_gb

    String bcftools_docker
  }

  String out_filename = out_prefix + ".vcf.gz"

  Int default_disk_gb = ceil(2.5 * size(vcfs, "GB")) + 10

  command <<<
    set -eu -o pipefail

    bcftools concat \
      ~{bcftools_concat_options} \
      --file-list ~{write_lines(vcfs)} \
      -O z \
      -o ~{out_filename} \
      --threads ~{cpu_cores}

    tabix -p vcf -f ~{out_filename}
  >>>

  output {
    File merged_vcf = "~{out_filename}"
    File merged_vcf_idx = "~{out_filename}.tbi"
  }

  runtime {
    docker: bcftools_docker
    memory: mem_gb + " GB"
    cpu: cpu_cores
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    preemptible: 3
  }
}


task GetSamplesFromVcfHeader {
  input {
    File vcf
    File vcf_idx
    String bcftools_docker
  }

  String out_filename = basename(vcf, ".vcf.gz") + ".samples.list"
  Int disk_gb = ceil(1.2 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    bcftools query -l ~{vcf} > ~{out_filename}
  >>>

  output {
    File sample_list = out_filename
    Int n_samples = length(read_lines(out_filename))
  }

  runtime {
    docker: bcftools_docker
    memory: "3.75 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


# Make an empty text file as a placeholder
task MakeDummyFile {
  input {}

  command <<<
    touch "dummy.txt"
  >>>

  output {
    File empty_file = "dummy.txt"
  }

  runtime {
    docker: "ubuntu:plucky-20251001"
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk 15 HDD"
    preemptible: 3
    maxRetries: 2
  }
}


# Make an empty benchmark bed to backfill for missing values in optional arrays
task MakeEmptyBenchBed {
  input {
    String docker
    String header_format = "site"
  }

  command <<<
    set -eu -o pipefail

    if [ "~{header_format}" == "site" ]; then
      echo -e "#chrom\tstart\tend\tvid\tclass\tsubclass\tsize\tac\taf\tfreq_het\tfreq_hom\thwe" \
      | bgzip -c > dummy.bed.gz
    fi
  >>>

  output {
    File empty_bed = "dummy.bed.gz"
  }

  runtime {
    docker: docker
    memory: "1.5 GB"
    cpu: 1
    disks: "local-disk 20 HDD"
    preemptible: 3
  }
}


task MakeHeaderFiller {
  input {}

  command <<<
    set -eu -o pipefail

    echo "##FILTER=<ID=MULTIALLELIC,Description=\"Multiallelic site\">" > header.supp.vcf
    echo "##FORMAT=<ID=RD_CN,Number=1,Type=Integer,Description=\"Predicted copy state\">" >> header.supp.vcf
    echo "##INFO=<ID=AC_Het,Number=A,Type=Integer,Description=\"Heterozygous allele counts\">" >> header.supp.vcf
    echo "##INFO=<ID=AC_Hom,Number=A,Type=Integer,Description=\"Homozygous allele counts\">" >> header.supp.vcf
    echo "##INFO=<ID=AC_Hemi,Number=A,Type=Integer,Description=\"Hemizygous allele counts\">" >> header.supp.vcf
    echo "##INFO=<ID=CN_NONREF_COUNT,Number=1,Type=Integer,Description=\"Nondip count.\">" >> header.supp.vcf
    echo "##INFO=<ID=CN_NONREF_FREQ,Number=1,Type=Float,Description=\"CNV frequency\">" >> header.supp.vcf
    echo "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">" >> header.supp.vcf
    echo "##INFO=<ID=HWE,Number=A,Type=Float,Description=\"HWE test\">" >> header.supp.vcf
    echo "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">" >> header.supp.vcf
    echo "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" >> header.supp.vcf
    
  >>>

  output {
    File supp_vcf_header = "header.supp.vcf"
  }

  runtime {
    docker: "ubuntu:plucky-20251001"
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk 20 HDD"
    preemptible: 3
    max_retries: 2
  }
}


# Duplicated from Utilities.wdl
task MakeTabixIndex {
  input {
    File input_file
    String file_type = "vcf"
    String docker
    Int n_cpu = 2
    Float mem_gb = 3.5
  }

  String outfile = basename(input_file, "gz") + "gz.tbi"
  Int disk_gb = ceil(1.25 * size(input_file, "GB")) + 5

  command <<<
    set -eu -o pipefail

    # Due to weirdness with call cacheing and where tabix writes 
    # indexes by default, it seems safest to relocate the input
    # file to the pwd before indexing it
    
    mv ~{input_file} ./
    tabix -p ~{file_type} -f ~{basename(input_file)}
  >>>

  output {
    File tbi = "~{outfile}"
  }

  runtime {
    docker: docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}


# Checks a VCF header for fields indicating that it might contain mCNVs
# Functionally equivalent to StreamedMcnvHeaderCheck, but localizes entire VCF
# Unable to disable localization_optional in parameter_meta so we must 
# have two versions of roughly the same function
task McnvHeaderCheck {
  input {
    File vcf
    File vcf_idx
    String docker
  }

  Int disk_gb = ceil(1.5 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    bcftools view --header-only ~{vcf} > header.vcf

    # Check for header rows potentially indicative of MCNVs
    touch mcnv.header.vcf
    fgrep "ALT=<ID=CNV" header.vcf > mcnv.header.vcf || true
    fgrep "ALT=<ID=MCNV" header.vcf >> mcnv.header.vcf || true
    fgrep "FILTER=<ID=MULTIALLELIC" header.vcf >> mcnv.header.vcf || true
    if [ $( cat mcnv.header.vcf | wc -l ) -gt 0 ]; then
      echo "true" > has_mcnvs.txt
    else
      echo "false" > has_mcnvs.txt
    fi
  >>>

  output {
    Boolean has_mcnvs = read_boolean("has_mcnvs.txt")
  }

  runtime {
    docker: docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
    maxRetries: 2
  }
}


# Parse a GATK-style interval_list and output as Array[[index, interval_string]] for scattering
# It appears WDL read_tsv does not allow coercion to Pair objects, so we output as Array[Array[String]]
task ParseIntervals {
  input {
    File intervals_list
    String docker
  }

  command <<<
    set -eu -o pipefail

    fgrep -v "#" ~{intervals_list} | fgrep -v "@" | sort -V \
    | awk -v OFS="\t" '{ print NR, $0 }' > intervals.clean.tsv
  >>>

  output {
    Array[Array[String]] interval_info = read_tsv("intervals.clean.tsv")
  }

  runtime {
    docker: docker
    memory: "2 GB"
    cpu: 1
    disks: "local-disk 25 HDD"
    preemptible: 3
  }
}


# Preps a sites.bed file for site comparison by generating:
# 1. A strict "query" set, which only includes variants with POS within an 
# eval interval, >50% coverage by all eval intervals, and size within [min_size, max_size]
# 2. A lenient "ref" set, which includes any variant overlapping any eval_interval
# and inclusive of all variants with size within [min_size / 3, 3 * max_size]
task PrepSites {
  input {
    Array[File] beds
    Array[File] bed_idxs
    File eval_interval_bed

    Int min_size = 0
    Int max_size = 1000000000
    Float lenient_size_scalar = 3.0
    Float strict_interval_coverage = 0.5
    
    String prefix
    
    String g2c_analysis_docker
  }

  Int loose_min_size = floor(min_size / lenient_size_scalar)
  Int loose_max_size = ceil(lenient_size_scalar * max_size)

  Int disk_gb = ceil(3 * size(beds, "GB")) + 10

  command <<<
    set -eu -o pipefail

    # Link bed indexes if necessary
    while read bed_path idx_path; do
      if [ "$idx_path" != "$bed_path.tbi" ]; then
        ln -s "$idx_path" "$bed_path.tbi"
      fi
    done < <( paste ~{write_lines(beds)} ~{write_lines(bed_idxs)} )

    # Save header line for first BED
    tabix --only-header -p bed -R ~{eval_interval_bed} ~{beds[0]} > header.bed

    # Loop over each bed and use tabix to extract only the variants needed
    while read bed_path; do
      tabix -p bed -R ~{eval_interval_bed} $bed_path
    done < ~{write_lines(beds)} \
    | awk -v FS="\t" -v min_size=~{loose_min_size} -v max_size=~{loose_max_size} \
      '{ if ($7>=min_size && $7<=max_size) print }' \
    | sort -Vk1,1 -k2,2n -k3,3n | uniq \
    | cat header.bed - | bgzip -c \
    > ~{prefix}.ref.bed.gz

    # Further filter the lenient "ref" set to produce a strict "query" set
    if [ $( zcat ~{prefix}.ref.bed.gz | fgrep -v "#" | wc -l ) -gt 0 ]; then
      /opt/pancan_germline_wgs/scripts/qc/vcf_qc/enforce_strict_intervals.py \
        -i ~{prefix}.ref.bed.gz \
        -t ~{eval_interval_bed} \
        -f ~{strict_interval_coverage} \
        -m ~{min_size} \
        -M ~{max_size} \
        -o ~{prefix}.query.bed.gz
    else
      cp ~{prefix}.ref.bed.gz ~{prefix}.query.bed.gz
    fi
  >>>

  output {
    File ref_bed = "~{prefix}.ref.bed.gz"
    File query_bed = "~{prefix}.query.bed.gz"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}


# Divide a BED-style intervals file into shards balanced on total genomic footprint
task ShardIntervals {
  input {
    File intervals_bed
    Int n_shards
    String prefix
    String g2c_analysis_docker

    Float mem_gb = 3.75
    Int n_cpu = 2
    Int? disk_gb
  }

  Int default_disk_gb = ceil(10 * size(intervals_bed, "GB")) + 5
  Int n_shards_use = if n_shards < 1 then 1 else n_shards

  command <<<
    set -eu -o pipefail

    # Ensure intervals file is simple, sorted, bgzipped, and merged
    zcat ~{intervals_bed} \
    | fgrep -v "#" \
    | cut -f1-3 \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | bedtools merge -i - \
    | bgzip -c \
    > intervals.prepped.bed.gz

    # Estimate total size per shard
    shard_size=$( zcat intervals.prepped.bed.gz \
                  | awk -v denom=~{n_shards_use} \
                    '{ sum+=$3-$2 }END{ printf "%.0f", sum / denom }' )

    # Ensure no interval is larger than shard_size
    /opt/pancan_germline_wgs/scripts/gatkhc_helpers/split_intervals.py \
      -i intervals.prepped.bed.gz \
      -t $shard_size \
      --max-interval-size $shard_size \
      --bed-style \
      -o intervals.prepped.split.bed
    bgzip intervals.prepped.split.bed

    # Shard intervals
    k=1
    bp_remain=$shard_size
    while read chrom start end; do

      size=$(( $end - $start ))

      if [ $size -le $bp_remain ]; then
        echo -e "$chrom\t$start\t$end" >> "~{prefix}.$k.bed"
        bp_remain=$(( $bp_remain - $size ))

      else
        new_end=$(( $start + $bp_remain ))
        echo -e "$chrom\t$start\t$new_end" >> "~{prefix}.$k.bed"
        
        # Because we split all intervals to be <= shard_size above,
        # we know that the remainder must always be < shard_size
        ((k++))
        start=$new_end
        size=$(( $end - $start ))
        echo -e "$chrom\t$start\t$end" >> "~{prefix}.$k.bed"
        bp_remain=$(( $shard_size - $size ))
      fi

    done < <( zcat intervals.prepped.split.bed.gz )
    
    find ./ -name "~{prefix}.*.bed" | xargs -I {} bgzip {}
  >>>

  output {
    Array[File] interval_shards = glob("~{prefix}.*.bed.gz")
  }

  runtime {
    docker: g2c_analysis_docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) +" HDD"
    preemptible: 1
    maxRetries: 1
  }
}


task ShardTextFile {
  input {
    File input_file
    Int n_splits
    String out_prefix
    Boolean shuffle = false
    Float mem_gb = 3.5
    String g2c_analysis_docker
  }

  Int disk_gb = ceil(3 * size(input_file, "GB")) + 10
  String shuffle_cmd = if shuffle then "--shuffle" else ""

  command <<<
    set -eu -o pipefail

    if [ ~{n_splits} -gt 1 ]; then
      /opt/pancan_germline_wgs/scripts/utilities/evenSplitter.R \
        -S ~{n_splits} \
        ~{shuffle_cmd} \
        ~{input_file} \
        ~{out_prefix}
    else
      cp ~{input_file} "~{out_prefix}1"
    fi
  >>>

  output {
    Array[File] shards = glob("~{out_prefix}*")
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "~{mem_gb} GB"
    cpu: 2
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}


# Duplicated from Utilities.wdl
task ShardVcf {
  input {
    File vcf
    File vcf_idx
    Int records_per_shard
    String bcftools_docker
    Int? disk_gb
    Int n_preemptible = 3
  }

  String out_prefix = basename(vcf, ".vcf.gz") + ".sharded"
  Int use_disk_gb = select_first([disk_gb, ceil(5 * size(vcf, "GB")) + 20])

  command <<<
    set -eu -o pipefail

    # Make an empty shard in case the input VCF is totally empty
    bcftools view -h ~{vcf} | bgzip -c > "~{out_prefix}.0.vcf.gz"

    bcftools +scatter \
      -O z3 -o . -p "~{out_prefix}". \
      -n ~{records_per_shard} \
      ~{vcf}

    # Print all VCFs to stdout for logging purposes
    find ./ -name "*.vcf.gz"

    # Index all shards
    find ./ -name "~{out_prefix}.*.vcf.gz" \
    | xargs -I {} tabix -p vcf -f {}
  >>>

  output {
    Array[File] vcf_shards = glob("~{out_prefix}.*.vcf.gz")
    Array[File] vcf_shard_idxs = glob("~{out_prefix}.*.vcf.gz.tbi")
  }

  runtime {
    cpu: 2
    memory: "3.75 GiB"
    disks: "local-disk " + use_disk_gb + " HDD"
    bootDiskSizeGb: 10
    docker: bcftools_docker
    preemptible: n_preemptible
    maxRetries: 1
  }
}


# Extract a slice from a VCF with remote indexing rather than localizing the entire file
task StreamSliceVcf {
  input {
    File vcf
    File vcf_idx
    String interval

    File? samples_list

    String? outfile_name
    
    Int? disk_gb
    Float mem_gb = 4
    Int n_cpu = 2

    String bcftools_docker    
  }

  String outfile = select_first([outfile_name, basename(vcf, ".vcf.gz") + ".sliced.vcf.gz"])

  String samples_cmd = if defined(samples_list) then "--samples-file ~{basename(select_first([samples_list]))} --force-samples" else ""

  Int default_disk_gb = ceil(size(vcf, "GB")) + 10
  Int hdd_gb = select_first([disk_gb, default_disk_gb])

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

    # Localize samples list to pwd if provided
    if [ ~{defined(samples_list)} ]; then
      cp ~{samples_list} ./
    fi

    # Stream VCF to interval of interest
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
    bcftools view --regions "~{interval}" ~{samples_cmd} ~{vcf} \
    | awk -v min_pos="$min_pos" -v max_pos="$max_pos" \
      '{ if ($1 ~ "^#" || ($2 >= min_pos && $2 <= max_pos)) print }' \
    | bcftools view -Oz -o "~{outfile}"

    # Index slice with tabix
    tabix -p vcf -f ~{outfile}
  >>>

  output {
    File vcf_slice = "~{outfile}"
    File vcf_slice_idx = "~{outfile}.tbi"
  }

  runtime {
    docker: bcftools_docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + hdd_gb + " HDD"
    preemptible: 3
  }
}


task SubsetVcfByVids {
  input {
    File vcf
    File vcf_idx
    File vids_list
    String bcftools_docker
    Int? disk_gb
    Int n_preemptible = 3
  }

  String out_fname = basename(vcf, ".vcf.gz") + ".subset.vcf.gz"
  Int use_disk_gb = select_first([disk_gb, ceil(2 * size(vcf, "GB")) + 10])

  command <<<
    set -eu -o pipefail

    bcftools view \
      --no-update \
      --include "ID=@~{vids_list}" \
      -Oz -o "~{out_fname}" \
      ~{vcf}
    tabix -p vcf -f "~{out_fname}"
  >>>

  output {
    File subsetted_vcf = out_fname
    File subsetted_vcf_idx = "~{out_fname}.tbi"
  }

  runtime {
    cpu: 1
    memory: "2 GiB"
    disks: "local-disk " + use_disk_gb + " HDD"
    bootDiskSizeGb: 10
    docker: bcftools_docker
    preemptible: n_preemptible
    maxRetries: 1
  }
}


# Sum one or more compressed distributions generated by CollectSiteMetrics
task SumCompressedDistribs {
  input {
    Array[File] distrib_tsvs
    Int n_key_columns = 2
    String out_prefix
    String g2c_analysis_docker
  }

  Int disk_gb = ceil(1.5 * size(distrib_tsvs, "GB")) + 10

  command <<<
    set -eu -o pipefail
    
    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/sum_compressed_distribs.py \
      -o "~{out_prefix}.merged.tsv" \
      -k ~{n_key_columns} \
      ~{sep=" " distrib_tsvs}
    gzip -f "~{out_prefix}.merged.tsv"
  >>>

  output {
    File merged_distrib = "~{out_prefix}.merged.tsv.gz"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}


