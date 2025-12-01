#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Fixes a prematurely truncated VCF by dropping the final truncated record
"""


import argparse
import gzip
import io
from g2cpy import determine_filetype
from sys import stdin, stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input-vcf', metavar='VCF', default='stdin',
                        help='Input .vcf. Stdin must be decompressed. Default: stdin')
    parser.add_argument('-o', '--output-vcf', metavar='VCF', default='stdout',
                        help='Path to output .vcf. Note that due to the nature ' +
                        'of this operation, this output cannot be automatically ' +
                        'compressed; thus, it is strongly encouraged to pipe to ' +
                        'bgzip. Default: stdout')
    args = parser.parse_args()

    # Open connection to output VCF
    if args.output_vcf in 'stdout /dev/stdout -'.split():
        outvcf = stdout
    else:
        outvcf = open(args.output_vcf, 'w')

    # Iterate over each line of input VCF and, if complete, write to output VCF
    if 'compressed' in determine_filetype(args.input_vcf):
        with open(args.input_vcf, 'rb') as fin:
            with gzip.GzipFile(fileobj=fin) as fgz:
                for _ in (True,):
                    with io.TextIOWrapper(fgz, encoding='utf-8') as ft:
                        try:
                            for line in ft:
                                if line.endswith('\n'):
                                    outvcf.write(line)
                                else:
                                    break
                        except:
                            break
    else:
        if args.input_vcf in 'stdin /dev/stdin -'.split():
            invcf = stdin
        else:
            invcf = open(args.input_vcf)
        for line in invcf.readlines():
            if line.endswith('\n'):
                outvcf.write(line)
            else:
                break

    # Close output VCF to flush buffer
    outvcf.close()

if __name__ == '__main__':
    main()

