#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Mask SV GTs based on quality prior to training an SNV imputation model
"""


import argparse
import numpy as np
import pysam
from g2cpy import is_multiallelic, parse_gt
from sys import stdin, stdout


def mask_gts(record, key='SL', always_keep_ref=True):
    """
    Mask genotypes with FORMAT/`key` < mean(FORMAT/`key`)
    """

    vals, sids = [], []
    for sid, sdat in record.samples.items():
        gtd = parse_gt(sdat.get('GT'))

        # Typically we don't need to worry about missing or ref GTs (unless optioned)
        if gtd['AN'] == 0 \
        or (always_keep_ref and gtd['AC'] == 0):
            continue

        # Otherwise, retain this sample's information
        vals.append(sdat.get(key, np.nan))
        sids.append(sid)
    
    # For the purposes of excluding low-quality genotypes, we mask anything with
    # quality score below the mean of all quality scores. This is intentionally not
    # robust to outliers, as we actually _care_ about outliers in this case.
    # To avoid penalizing underdispersed integer quality score distributions, we
    # check if all quality values are integers. If so, we take the floor of `cutoff`
    cutoff = np.nanmean(vals)
    def _is_int_like(v):
        if v is None or np.isnan(v):
            return True
        try:
            int(v)
        except:
            return False
        return True
    if all([_is_int_like(v) for v in vals]):
        cutoff = int(np.floor(cutoff))
    for sid in sids:
        v = record.samples[sid].get(key, None)
        if v is not None:
            if v < cutoff:
                gtd = parse_gt(record.samples[sid]['GT'])
                nocall = tuple([None for k in range(gtd['AN'] + gtd['N_missing'])])
                record.samples[sid]['GT'] = nocall

    return record


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input-vcf', default='stdin', metavar='VCF', 
                        type=str, help='Input .vcf. [default: stdin].')
    parser.add_argument('-o', '--output-vcf', default='stdout', metavar='VCF',
                        type=str, help='Output .vcf. [default: stdout]')
    parser.add_argument('-q', '--quality-field', type=str, default='SL',
                        help='FORMAT field to use for automated quality filtering')
    args = parser.parse_args()

    # Open connection to input VCF
    if args.input_vcf in 'stdin /dev/stdin -'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.input_vcf)

    # Open connection to output VCF
    if args.output_vcf in 'stdout /dev/stdout -'.split():
        outvcf = pysam.VariantFile(stdout, 'w', header=invcf.header.copy())
    else:
        outvcf = pysam.VariantFile(args.output_vcf, 'w', header=invcf.header.copy())

    # Iterate over input VCF and dynamically mask each record
    for record in invcf:

        # Mask genotypes for biallelic variants
        if not is_multiallelic(record):
            record = mask_gts(record, args.quality_field)

        # Write to output VCF
        outvcf.write(record)

    # Close connection to output VCF
    outvcf.close()


if __name__ == '__main__':
    main()

