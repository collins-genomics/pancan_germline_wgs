#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Merge HMF data into a set of Terra-style manifests
"""


import argparse
import csv
import json
import pandas as pd


def load_superclass_table(supermap_tsv):
    """
    Build a dictionary for classifying tumors into G2C cancer types
    """

    supermap = {}

    superdf = pd.read_csv(supermap_tsv, sep='\t')

    for idx, rvals in superdf.iterrows():

        pl = str(rvals.primaryTumorLocation)
        sl = str(rvals.primaryTumorSubLocation)
        pt = str(rvals.primaryTumorType)
        g2c_class = rvals['class']

        if pl not in supermap.keys():
            supermap[pl] = {}

        if sl not in supermap[pl].keys():
            supermap[pl][sl] = {}

        if pt not in supermap[pl][sl].keys():
            supermap[pl][sl][pt] = g2c_class

    return supermap


def map_g2c_cancer_type(vals, supermap):
    """
    Converts a single HMF metadata entry to a G2C cancer type
    """

    pl = str(vals.primaryTumorLocation)
    sl = str(vals.primaryTumorSubLocation)
    pt = str(vals.primaryTumorType)

    try:
        return supermap[pl][sl][pt]
    except:
        import pdb; pdb.set_trace()


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--metadata', help='HMF metadata.tsv',
                        required=True)
    parser.add_argument('--manifest', help='HMF manifest.json', 
                        required=True)
    parser.add_argument('--tumor-superclass-table', help='.tsv with instructions ' +
                        'for classifying tumors based on primaryTumorLocation, ' +
                        'primaryTumorSubLocation, and primaryTumorType.')
    parser.add_argument('--join', default='outer', help='behavior for joining ' +
                        '--metadata and --manifest.')
    parser.add_argument('-o', '--outdir', help='path to output directory',
                        default='./')
    args = parser.parse_args()

    # Load metadata (do not clean yet)
    md = pd.read_csv(args.metadata, sep='\t')

    # Annotate metadata with G2C cancer classification, if optioned
    supermap = {}
    if args.tumor_superclass_table is not None:
        supermap = load_superclass_table(args.tumor_superclass_table)
        md.loc[:, 'g2c_cancer_type'] = md.apply(map_g2c_cancer_type, axis=1, supermap=supermap)
    
    # Parse file manifest .json and convert all entries into pd.DataFrame
    mfst_cols = 'sampleId tumor_cram tumor_crai normal_cram normal_crai ' + \
                'germline_vcf germline_vcf_idx has_rnaseq'
    mfst = pd.DataFrame(columns = mfst_cols.split())
    with open(args.manifest) as fin:
        for sdat in json.load(fin)['data']:
            sid = sdat.get('sampleId')
            tcram = sdat.get('crams', []).get('tumor', []).get('url')
            tcrai = sdat.get('crams', []).get('tumorIndex', []).get('url')
            ncram = sdat.get('crams', []).get('normal', []).get('url')
            ncrai = sdat.get('crams', []).get('normalIndex', []).get('url')
            vcf = sdat.get('germline', []).get('variants', []).get('url')
            vcf_idx = sdat.get('germline', []).get('index', []).get('url')
            has_rna = len(sdat.get('rna', []).get('fastq', [])) > 0
            srow = [sid, tcram, tcrai, ncram, ncrai, vcf, vcf_idx, has_rna]
            mfst.loc[len(mfst.index)] = srow

    # Merge metadata and manifest
    df = md.merge(mfst, how=args.join, on='sampleId')

    # Write full metadata as tsv to outdir
    df.to_csv(args.outdir + '/HMF.combined_metadata_and_manifest.tsv.gz',
              sep='\t', index=False, na_rep='NA')

    # Sort based on tumor biopsy date (oldest first) before deduplicating on
    # patient ID such that only one T/N pair is retained per patient
    df = df.sort_values('biopsyDate').\
            drop_duplicates(subset='hmfPatientId', keep='first')
    df.to_csv(args.outdir + '/HMF.combined_metadata_and_manifest.unique_patients.tsv.gz',
              sep='\t', index=False, na_rep='NA')

    # Write terra-style metadata for WGS processing to outdir
    twgs_keeper_cols = 'hmfPatientId primaryTumorLocation primaryTumorSubLocation ' + \
                       'primaryTumorType normal_cram normal_crai'
    if args.tumor_superclass_table is not None:
        twgs_keeper_cols += ' g2c_cancer_type'
        keep_rows = (~df.normal_cram.isna()) & (df.g2c_cancer_type != 'Other')
    else:
        keep_rows = ~df.normal_cram.isna()
    twgs_keeper_cols = list(set(twgs_keeper_cols.split()).intersection(set(df.columns)))
    twgs = df.loc[keep_rows, twgs_keeper_cols]
    twgs.rename(columns = {'hmfPatientId' : 'entity:sample_id',
                           'normal_cram' : 'hg19_cram',
                           'normal_crai' : 'hg19_crai'},
                inplace=True)
    twgs.to_csv(args.outdir + '/HMF.wgs_processing.terra_manifest.tsv',
                sep='\t', index=False, na_rep='NA')

    # Write terra-style metadata for hg19 germline SNV/indel processing to outdir
    tvcf_keeper_cols = 'hmfPatientId primaryTumorLocation primaryTumorSubLocation ' + \
                       'primaryTumorType germline_vcf germline_vcf_idx'
    if args.tumor_superclass_table is not None:
        tvcf_keeper_cols += ' g2c_cancer_type'
        keep_cancers = 'Lung Colorectal Pancreas Kidney Prostate Breast'.split()
        keep_rows = (~df.germline_vcf.isna()) & (df.g2c_cancer_type.isin(keep_cancers))
    else:
        keep_rows = ~df.germline_vcf.isna()
    tvcf_keeper_cols = list(set(tvcf_keeper_cols.split()).intersection(set(df.columns)))
    tvcf = df.loc[keep_rows, tvcf_keeper_cols]
    tvcf.rename(columns = {'hmfPatientId' : 'entity:sample_id'},
                inplace=True)
    tvcf.to_csv(args.outdir + '/HMF.hg19_germline_vcf_processing.terra_manifest.tsv',
                sep='\t', index=False, na_rep='NA')


if __name__ == '__main__':
    main()

