#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Reassign GQ as a transformation of SL
"""


import argparse
import numpy as np
import pysam
from g2cpy import parse_gt
from sys import stdin, stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--in-vcf', metavar='vcf', default='stdin',
                        help='input .vcf [default: stdin]')
    parser.add_argument('-o', '--out-vcf', metavar='vcf', default='stdout',
                        help='output .vcf [default: stdout]')
    args = parser.parse_args()

    # Open connection to input vcf
    if args.in_vcf in 'stdin /dev/stdin -'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.in_vcf)

    # Get list of samples
    samples = [s for s in invcf.header.samples]

    # Open connection to output vcf
    if args.out_vcf in 'stdout /dev/stdout -'.split():
        outvcf = pysam.VariantFile(stdout, mode='w', header=invcf.header)
    else:
        outvcf = pysam.VariantFile(args.out_vcf, mode='w', header=invcf.header)

    # Process each record
    for record in invcf:

        gqs = []

        # Process each sample per record
        for sid in samples:
            sdat = record.samples[sid]

            # Skip no-call genotypes
            gt = parse_gt(sdat.get('GT', (None, None, )))
            if gt['N_missing'] == gt['AN']:
                continue

            # Only reassign GQ if SL is present and non-null
            if 'SL' not in sdat.keys():
                continue
            sl = sdat.get('SL')
            if sl is None or sl == '.':
                continue

            # SL reflects the odds that a genotype is non-reference
            # However, Phred reflects the probability of an *error*
            # Thus, we need to invert SL for all non-reference genotypes
            nonref = False
            if gt['AC'] > 0:
                nonref = True
                sl = -sl

            # Converting SL to GQ following the specifications on page 36 of the
            # All of Us Genomic Data Quality Report v7:
            # https://support.researchallofus.org/hc/en-us/articles/4617899955092-All-of-Us-Genomic-Quality-Report-ARCHIVED-C2022Q4R9-CDR-v7
            p = 1 / (1 + ((0.52/0.48) ** -sl))
            raw_gq = -10 * np.log10(p)
            gq = int(np.round(np.nanmax([np.nanmin([raw_gq, 99]), 0]), 0))

            # Assign GQ
            record.samples[sid]['GQ'] = gq
            if nonref:
                gqs.append(gq)

        # If any nonref GQs were recalibrated, reassign QUAL as 10 * mean nonref GQ
        if len(gqs) > 0:
            qual = 10 * np.mean(gqs)
            record.qual = int(np.round(qual, 0))

        # Write updated record to output VCF
        outvcf.write(record)

    # Clear buffers
    outvcf.close()


if __name__ == '__main__':
    main()


