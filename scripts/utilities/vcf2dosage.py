#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Convert a VCF into a simple allele dosage matrix
"""


import argparse
import pysam
from g2cpy import name_record, parse_gt
from sys import stdin, stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vcf-in', help='input .vcf', metavar='vcf', default='stdin')
    parser.add_argument('--tsv-out', help='output .tsv', metavar='.tsv', default='stdout')
    args = parser.parse_args()

    # Open connection to input vcf
    if args.vcf_in in '- stdin /dev/stdin'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.vcf_in)
    samples = [s for s in invcf.header.samples]

    # Open connection to output file
    if args.tsv_out in '- stdout /dev/stdout':
        outfile = stdout
    else:
        outfile = open(args.tsv_out, 'w')
    outfile.write('\t'.join(['#VID'] + samples) + '\n')

    # Iterate over records in invcf, convert to dosage vectors, and write to outfile
    for record in invcf.fetch():
        vid = record.id
        if vid is None:
            vid = name_record(record)
        outvals = [vid]
        for sid, sgt in record.samples.items():
            try:
                ac, an, nmiss = parse_gt(sgt['GT']).values()
                if an == 0:
                    dos = 'NA'
                else:
                    dos = str(ac)
            except:
                dos = 'NA'
            outvals.append(dos)
        outfile.write('\t'.join(outvals) + '\n')

    # Close connection to output file to clear buffer
    outfile.close()


if __name__ == '__main__':
    main()

