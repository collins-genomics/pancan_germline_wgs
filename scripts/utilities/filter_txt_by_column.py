#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Filter a text file based on intersecting one column with a list of eligible keys
"""


import argparse
import gzip
from sys import stdin, stdout
from g2cpy import determine_filetype


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('infile', help='Input file. Can be gzip compressed. ' +
                        'Also accepts stdin, /dev/stdin, and -')
    parser.add_argument('keyfile', help='List of eligible keys to include')
    parser.add_argument('-c', '--column-number', default=1, type=int,
                        help='Column number from infile to filter by [default: 1]')
    parser.add_argument('-d', '--delimiter', default='\t', type=str,
                        help='Delimiter in infile [default: tab]')
    parser.add_argument('-o', '--outfile', default='stdout',
                        help='Path to output file [default: stdout]')
    args = parser.parse_args()

    # Load keys into memory
    with open(args.keyfile) as fin:
        keys = set([k.rstrip() for k in fin.readlines()])

    # Open connection to input file
    if args.infile in 'stdin /dev/stdin -'.split():
        infile = stdin
    elif 'compressed' in determine_filetype(args.infile):
        infile = gzip.open(args.infile, 'rt')
    else:
        infile = open(args.infile)

    # Open connection to output file
    if args.outfile in 'stdout /dev/stdout -'.split():
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Filter input file
    for line in infile.readlines():
        # Always pass header line
        if line.startswith('#'):
            outfile.write(line)
        # Otherwise, split by delimiter and query vs. keys
        if line.split(args.delimiter)[args.column_number - 1] in keys:
            outfile.write(line)

    # Close file handles to clear buffer
    infile.close()
    outfile.close()


if __name__ == '__main__':
    main()
