#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Helper script for Terra/Terra-like implementations of the G2C VCF QC workflow
"""


import json
from sys import argv


def main():
    infile = argv[1]
    outfile = argv[2]

    # Read contents of input array
    with open(infile) as fin:
        jdat = json.load(fin)

    # Determine dimensions, assuming all axes are rectangular (same dimensions across strata)
    n_shards = len(jdat)
    try:
        n_datasets = len(jdat[0])
    except:
        n_datasets = 0
    try:
        n_intervals = len(jdat[0][0])
    except:
        n_intervals = 0

    # Populate new nested list with correct structure for QC pipeline
    jout = []
    for i in range(n_datasets):
        jout.append([])
        for j in range(n_intervals):
            jout[i].append([])
            for k in range(n_shards):
                jout[i][j].append(jdat[k][i][j])

    # Write transposed array to output file
    with open(outfile, 'w') as fout:
        json.dump(jout, fout)


if __name__ == "__main__":
    main()
