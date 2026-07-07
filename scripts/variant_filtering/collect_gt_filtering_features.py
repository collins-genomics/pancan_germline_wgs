#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Gather genotype-level annotations for GT filtering model
Input VCF must strictly contain only a single variant type (snv, indel, or SVs)
"""


import argparse
import gzip
import g2cpy
import pysam
from sys import stdin, stdout


def gather_gt_features(sdat):
    """
    Main function to collect genotype-level features
    """

    outvals, outcols = [], []

    # Collect all numeric features
    # TODO: implement this

    # Tokenize all non-numeric features
    # TODO: implement this

    return outvals, outcols


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input-vcf', default='stdin', metavar='VCF', 
                        type=str, help='Input .vcf')
    parser.add_argument('-o', '--output-tsv', default='stdout', metavar='TSV',
                        type=str, help='Output .tsv')
    parser.add_argument('-a', '--collect-all', action='store_true',
                        help='Collect all genotypes (default: only collect non-ref)')
    parser.add_argument('-p', '--precision', default=3, metavar='integer',
                        type=int, help='Floating point precision')
    args = parser.parse_args()

    # Open connection to input VCF
    if args.input_vcf in 'stdin /dev/stdin -'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.input_vcf)

    # Open connection to output TSV
    if args.output_tsv in 'stdin /dev/stdin -'.split():
        outfile = stdout
    else:
        if 'compressed' in g2cpy.determine_filetype(args.output_tsv):
            outfile = gzip.open(args.output_tsv, 'wt', encoding='utf-8')
        else:
            outfile = open(args.output_tsv, 'w')

    # Iterate over input VCF records and samples within each record
    header_written = False
    for i, record in enumerate(invcf):
        for sid, sdat in record.samples.items():

            # Check if this sample needs to be processed
            sgt = g2cpy.parse_gt(sdat['GT'])
            if sgt['AC'] == 0 and not args.collect_all:
                continue

            # Collect basic info
            vid = record.id
            outvals = [vid, sid, sgt['AC'], sgt['AN'] - sgt['N_missing']]
            if i == 0:
                outcols = 'vid sid sac san'.split()

            # Collect other features
            add_outvals, add_outcols = gather_gt_features(sdat)
            outvals += add_outvals
            outcols += add_outcols

            # Write header if first line
            if i == 0 and not header_written:
                outfile.write('\t'.join(outcols) + '\n')
                header_written = True

            # Format & write output line
            outline = []
            for v in outvals:
                if isinstance(v, float):
                    outline.append(f"{v:.{args.precision}e}")
                else:
                    outline.append(str(v))
            outfile.write('\t'.join(outline) + '\n')

    # Close connection to output TSV
    outfile.close()


if __name__ == '__main__':
    main()

