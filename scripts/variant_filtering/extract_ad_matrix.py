#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Extract allele dosages from a GATK-formatted VCF
Returns a matrix of samples (rows) by variants (columns)
"""


import argparse
import pandas as pd
import pysam
from g2cpy import parse_gt, name_record


def get_ad(sdat, cov=None, ploidy=2):
    """
    Extracts allele dosage from a GATK-style FORMAT entry
    """

    # If GT is missing, return None
    if all([a is None for a in sdat.get('GT', (None, ))]):
        return None

    # Wherever possible, use quantitative estimate of allele dosage
    ad = sdat.get('AD')
    dp = sdat.get('DP')
    dos = None
    if ad is not None:
        nonref_ad = sum(ad[1:])
        if cov is not None:
            if dp is not None:
                ad_denom = dp / cov
            else:
                ad_denom = cov
        else:
            ad_denom = dp
        if ad_denom is not None:
            if ad_denom > 0:
                dos = ploidy * nonref_ad / ad_denom
    
    # Otherwise, just use simple GT field parsing, which should ~always be defined
    if dos is None:
        gt_parts = parse_gt(sdat.get('GT', (None, None, )))
        if gt_parts['N_missing'] == gt_parts['AN']:
            return None
        else:
            dos = gt_parts['AC']

    return float(dos)


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input-vcf', required=True, metavar='VCF', 
                        help='Input .vcf. Required.')
    parser.add_argument('-v', '--variant-ids', metavar='.txt', 
                        help='Optional list of variant IDs to retain')
    parser.add_argument('-c', '--coverage-tsv', metavar='.tsv', 
                        help='Optional two-column .tsv of sample IDs and median ' +
                             'diploid coverage values (for more accurate dosages).')
    parser.add_argument('-o', '--output-tsv', metavar='.tsv', default='stdout',
                        help='Path to output .tsv')
    args = parser.parse_args()

    # Open connection to input VCF
    invcf = pysam.VariantFile(args.input_vcf)

    # Load list of variant IDs to retain, if provided
    vids = None
    if args.variant_ids is not None:
        with open(args.variant_ids) as fin:
            vids = set([l.rstrip() for l in fin.readlines()])

    # Load a map of samples : median diploid coverage, if provided
    s_cov = dict()
    if args.coverage_tsv is not None:
        with open(args.coverage_tsv) as fin:
            s_cov = {s : float(c) for s, c in csv.reader(fin, delimiter='\t')}

    # Iterate over input VCF and route each SV to the correct output VCF
    ad_res = dict()
    for record in invcf:

        # Enforce --variant-ids if optioned
        if record.id is None:
            record.id = name_record(record)
        if vids is not None:
            if record.id not in vids:
                continue

        # Traverse over samples and compute allele dosages
        ad_map = dict()
        ad_res[record.id] = {sid : get_ad(sdat, s_cov.get(sid)) 
                             for sid, sdat in record.samples.items()}

    # Structure output as a matrix and write to --output-tsv
    out_df = pd.DataFrame.from_dict(ad_res)
    out_df.index.name = 'sample'
    out_df.to_csv(args.output_tsv, sep='\t', index=True, na_rep='NA')
    

if __name__ == '__main__':
    main()

