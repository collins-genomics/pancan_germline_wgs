#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Split a GATK-HC VCF into SNVs, small indels, and large indels
"""


import argparse
import numpy as np
import pysam
from g2cpy import classify_record
from sys import stdin, stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input-vcf', default='stdin', metavar='VCF', 
                        type=str, help='Input .vcf. [default: stdin].')
    parser.add_argument('-o', '--output-prefix', 
                        help='Prefix for output .vcfs [default: gatkhc]',
                        default='gatkhc', metavar='string', type=str)
    parser.add_argument('--large-indel-min-size', default=10, type=int, metavar='Int',
                        help='Minimum size for "large" indels [default: 10]')
    parser.add_argument('--largest-indel-log', default='stdout', type=str, 
                        metavar='txt', help='Output file for reporting the ' +
                        'length of the longest indel [default: stdout]')
    parser.add_argument('--make-large-indel-bed', action='store_true',
                        help='Also write out a sites .BED file for large indels')
    args = parser.parse_args()

    # Open connection to input VCF
    if args.input_vcf in 'stdin /dev/stdin -'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.input_vcf)
    
    # Open connection to output VCFs
    header = invcf.header.copy()
    out_vcf_fnames = {vc : '{}.{}.vcf.gz'.format(args.output_prefix, vc) 
                      for vc in 'snv small_indel large_indel'.split()}
    out_vcfs = {vc : pysam.VariantFile(fn, 'w', header=header) 
                for vc, fn in out_vcf_fnames.items()}

    # Open connection to output large indel sites BED, if optioned
    if args.make_large_indel_bed:
        li_bed_out = open('{}.large_indel.sites.bed'.format(args.output_prefix), 'w')

    # Iterate over input VCF and route each SV to the correct output VCF
    longest = 0
    for record in invcf:

        # Determine record routing
        vc, vsc, varlen = classify_record(record, return_varlen=True)
        if vc != 'snv':
            if np.abs(varlen) >= args.large_indel_min_size:
                vc = 'large_indel'
            else:
                vc = 'small_indel'
            
            # Update longest indel
            varlen = np.abs(varlen)
            if varlen > longest:
                longest = varlen

        # Write record to appropriate VCF
        out_vcfs[vc].write(record)

        # If optioned, write large indels to sites BED
        if vc == 'large_indel' and args.make_large_indel_bed:
            outline = '{}\t{}\t{}\n'.format(record.chrom, record.pos, record.pos+varlen)
            li_bed_out.write(outline)

    # Close connections to output VCFs
    for fh in out_vcfs.values():
        fh.close()

    # Close connection to large indel sites BED, if necessary
    if args.make_large_indel_bed:
        li_bed_out.close()

    # Report longest indel
    longest_line = '{}\n'.format(longest)
    if args.largest_indel_log in 'stdout /dev/stdout -'.split():
        stdout.write(longest_line)
    else:
        with open(args.largest_indel_log, 'w') as log_out:
            log_out.write(longest_line)


if __name__ == '__main__':
    main()

