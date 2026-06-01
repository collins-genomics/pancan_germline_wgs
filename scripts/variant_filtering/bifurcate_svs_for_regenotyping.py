#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Bifurcate an SV VCF to prepare for SNV-based GT refinement
"""


import argparse
import numpy as np
import pysam
from g2cpy import is_multiallelic, compute_allele_freq_stats
from sys import stdin, stdout


def label_sv(record, min_af=0, max_af=1, min_ac=0, max_ac=10e10, min_an=0):
    """
    Assign an SV to either 'eligible' or 'passthrough' VCF subsets based on frequency
    """

    # Only biallelic variants are eligible for regenotyping
    if is_multiallelic(record):
        return 'passthrough'

    # Get AF, AC, and AN if annotated, otherwise compute on the fly
    if all([f in record.info.keys() for f in 'AC AN AF'.split()]):
        ac = record.info.get('AC')[0]
        af = record.info.get('AF')[0]
        an = record.info.get('AN')
    else:
        freq_dat = compute_allele_freq_stats(record)
        ac = freq_dat.get('AC')
        af = freq.dat.get('AF')
        an = freq.dat.get('AN')

    # Compare frequency statistics to gates and emit corresponding label
    if af < min_af or af > max_af or ac < min_ac or ac > max_ac or an < min_an:
        return 'passthrough'
    else:
        return 'eligible'


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input-vcf', default='stdin', metavar='VCF', 
                        type=str, help='Input .vcf. [default: stdin].')
    parser.add_argument('-e', '--eligible-output-vcf', 
                        help='Output .vcf for eligible SVs. [default: stdout]',
                        default='stdout', metavar='VCF', type=str)
    parser.add_argument('-p', '--passthrough-output-vcf', 
                        help='Optional output .vcf for ineligible "pass-through" SVs',
                        metavar='VCF', type=str)
    parser.add_argument('--min-af', default=0.05, type=float, metavar='Float',
                        help='Lower AF threshold for eligibility [default: 0.05]')
    parser.add_argument('--max-af', default=0.95, type=float, metavar='Float',
                        help='Upper AF threshold for eligibility [default: 0.95]')
    parser.add_argument('--min-ac', default=20, type=int, metavar='Int',
                        help='Minimum AC for eligibility [default: 20]')
    parser.add_argument('--max-ac', default=10e10, type=int, metavar='Int',
                        help='Minimum AC for eligibility [default: ~infinite]')
    parser.add_argument('--min-an', default=100, type=int, metavar='Int',
                        help='Minimum AN for eligibility [default: 100]')
    args = parser.parse_args()

    # Open connection to input VCF
    if args.input_vcf in 'stdin /dev/stdin -'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.input_vcf)
    
    # Open connection to output VCF(s)
    header = invcf.header.copy()
    if args.eligible_output_vcf in 'stdout /dev/stdout -'.split():
        elig_vcf = pysam.VariantFile(stdout, 'w', header=header)
    else:
        elig_vcf = pysam.VariantFile(args.eligible_output_vcf, 'w', header=header)
    if args.passthrough_output_vcf is not None:
        pass_vcf = pysam.VariantFile(args.passthrough_output_vcf, 'w', header=header)

    # Iterate over input VCF and route each SV to the correct output VCF
    for record in invcf:

        # Determine record eligibility
        route = label_sv(record, args.min_af, args.max_af, args.min_ac, 
                         args.max_ac, args.min_an)
        if route == 'eligible':
            elig_vcf.write(record)
        elif args.passthrough_output_vcf is not None:
            pass_vcf.write(record)

    # Close connections to output VCF(s)
    elig_vcf.close()
    if args.passthrough_output_vcf is not None:
        pass_vcf.close()

if __name__ == '__main__':
    main()

