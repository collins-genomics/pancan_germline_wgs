#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Hail code to extract a desired subset of srWGS SNVs/indels from All of Us
Intended to be executed on a Dataproc cluster deployed in the Aou RW

Roughly follows the tutorial documentation available here:
https://workbench.researchallofus.org/workspaces/aou-rw-a5b0235e/howtoworkwithallofusgenomicdatahailplinkv8/analysis/preview/03_Manipulate%20Hail%20VariantDataset.ipynb
"""

# Import libraries
import hail as hl
from datetime import datetime
from sys import argv

# Set globals
main_bucket = 'gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45'
start = datetime.now()
vds_srwgs_path = 'gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/vds/hail.vds'
flagged_samples_path = 'gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/qc/flagged_samples.tsv'
keep_samples_path = main_bucket + '/refs/aou/AoU.G2C_samples_with_srWGS_sv.aou_ids.list'
example_header_url = 'gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/snv_indel/1KGP.srWGS.snv_indel.cleaned.chr1.vcf.gz'

# Read contig as only command-line positional argument
contig = argv[1]

# Initialize hail
print('\n\n\n\n\nBeginning SNV/indel curation for ' + contig + '\n')
hl.default_reference(new_default_reference='GRCh38')

# Read CDR v8 WGS .vds
print('Now loading VDS...')
vds = hl.vds.read_vds(vds_srwgs_path)

# Filter .vds to chromosome of interest
vds = hl.vds.filter_chromosomes(vds, keep=[contig])

# Filter .vds to non-flagged samples of interest
flagged_samples = hl.import_table(flagged_samples_path, key='s')
target_samples = hl.import_table(keep_samples_path, no_header=True).rename({'f0': 's'})
keep_samples_set = set(target_samples.s.collect()).difference(set(flagged_samples.s.collect()))
keep_samples = hl.Table.parallelize(
    [hl.struct(s=sample_id) for sample_id in keep_samples_set],
    key='s'
)
vds = hl.vds.filter_samples(vds, keep_samples, keep=True, remove_dead_alleles=True)

# Report number of variants and samples in initial subsetted VDS
msg = 'Loaded {:,} variants from {:,} samples for {} before filtering\n'
print(msg.format(*vds.variant_data.count(), contig))

# Convert VDS to dense MT
print('Now converting variants from VDS to dense MT...')
mt = vds.variant_data.annotate_entries(AD = hl.vds.local_to_global(vds.variant_data.LAD, 
                                                                   vds.variant_data.LA, 
                                                                   n_alleles=hl.len(vds.variant_data.alleles), 
                                                                   fill_value=0, number='R'))
mt = mt.annotate_entries(GT = hl.vds.lgt_to_gt(mt.LGT, mt.LA))
mt = mt.transmute_entries(FT = hl.if_else(mt.FT, "PASS", "FAIL"))
mt = hl.vds.to_dense_mt(hl.vds.VariantDataset(vds.reference_data, mt))
mt = mt.annotate_rows(info = hl.agg.call_stats(mt.GT, mt.alleles))

# Retain only FILTER PASS variants
print('Now filtering variants and genotypes...')
mt = mt.filter_rows(hl.len(mt.filters) == 0)

# Retain only FT PASS genotypes
mt = mt.filter_entries((mt.FT == "PASS") | hl.is_missing(mt.FT))

# Split mulitallelics
mt = hl.split_multi_hts(mt)

# Unphase genotypes
def unphasing_mt(mt):
    mt = mt.select_entries(GT = hl.unphased_diploid_gt_index_call(mt.GT.n_alt_alleles()))
    return mt
mt = unphasing_mt(mt)

# Drop all unnecessary annotations and update ones we want to keep
fields_to_drop_list = ['as_vets','as_vqsr', 'LAD', 'LGT', 'LA', 'tranche_data', 
                       'truth_sensitivity_snp_threshold', 
                       'truth_sensitivity_indel_threshold',
                       'snp_vqslod_threshold', 'indel_vqslod_threshold',
                       'was_split', 'GQ', 'PS', 'RGQ', 'FT', 'AD']
mt = mt.drop(*(f for f in fields_to_drop_list if f in mt.entry or f in mt.row or f in mt.col or f in mt.globals))
mt = mt.annotate_rows(
    info=mt.info.drop('homozygote_count')
)
mt = hl.variant_qc(mt)
mt = mt.annotate_rows(info=mt.info.annotate(AC=mt.variant_qc.AC[1:],
                                            AN=mt.variant_qc.AN,
                                            AF=mt.variant_qc.AF[1:]))

# Only retain variants where non-ref AC > 0 after all curation steps
mt = mt.filter_rows(hl.any(lambda x: x > 0, mt.variant_qc.AF[1:]))

# Not running the log below because it is quite slow; can be uncommented for debugging
# # Report number of variants and samples in output MT after filtering
# msg = 'Retained {:,} biallelic variants from {:,} samples for {} after filtering'
# print(msg.format(*mt.count(), contig))

# Convert MT to VCF & write to bucket
print('Now writing filtered variants as VCF...')
metadata = hl.get_vcf_metadata(example_header_url)
out_vcf = f'{main_bucket}/refs/aou/dense_vcfs/srwgs/snv_indel/AoU.srWGS.snv_indel.cleaned.{contig}.vcf.bgz'
hl.export_vcf(mt, out_vcf, tabix = False, metadata=metadata)

# Report completion
msg = 'Finished {} after {}'
stop = datetime.now()
print(msg.format(contig, str(stop-start)))
