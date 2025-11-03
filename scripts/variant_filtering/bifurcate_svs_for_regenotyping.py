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
import pysam
from g2cpy import is_multiallelic


def label_sv(record, min_af=0, max_af=1, min_ac=0, max_ac=10e10):
    """
    Assign an SV to either 'eligible' or 'passthrough' VCF subsets based on frequency
    """

    # Only biallelic variants are eligible for regenotyping
    if is_multiallelic(record):
        return 'passthrough'

    # Get AF and AC if annotated, otherwise compute on the fly
    af = record.info.get('AF', )
    ac = record.info.get('AC', )

    # Compare frequency statistics to gates and emit corresponding label
    if af < min_af or af > max_af:
        return 'passthrough'
    if ac < min_ac or ac > max_ac:
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
    parser.add_argument('-i', '--input-vcf', required=True, metavar='VCF', 
                        type=str, help='Input .vcf. Required.')
    parser.add_argument('-e', '--eligible-output-vcf', 
                        help='Output .vcf for eligible SVs',
                        required=True, metavar='VCF', type=str)
    parser.add_argument('-p', '--passthrough-output-vcf', 
                        help='Output .vcf for ineligible "pass-through" SVs',
                        metavar='VCF', type=str)
    parser.add_argument('--min-af', default=0.05, type=float, metavar='Float',
                        help='Lower AF threshold for eligibility [default: 0.05]')
    parser.add_argument('--max-af', default=0.95, type=float, metavar='Float',
                        help='Upper AF threshold for eligibility [default: 0.95]')
    parser.add_argument('--min-ac', default=20, type=int, metavar='Int',
                        help='Minimum AC for eligibility [default: 20]')
    parser.add_argument('--max-ac', default=10e10, type=int, metavar='Int',
                        help='Minimum AC for eligibility [default: ~infinite]')
    parser.add_argument('--samples', type=str, metavar='.txt',
                        help='List of sample IDs to include [default: keep all samples]')
    args = parser.parse_args()

    # Open connection to input VCF
    invcf = pysam.VariantFile(args.input_vcf)
    all_samples = set([s for s in invcf.header.samples])

    # Modify output header to only include --samples, if optioned
    header = invcf.header
    if args.samples is not None:
        # TODO: figure out how to do this
        print('DEV NOTE: sample subsetting not enabled yet')

    # Open connection to output VCF(s)
    elig_vcf = pysam.VariantFile(args.eligible_output_vcf, 'w', header=header)
    if args.passthrough_output_vcf is not None:
        pass_vcf = pysam.VariantFile(args.passthrough_output_vcf, 'w', header=header)

    # Iterate over input VCF and route each SV to the correct output VCF
    for record in invcf:

        # Determine record eligibility
        route = label_sv(record, args.min_af, args.max_af, args.min_ac, args.max_ac)

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

