#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Reshards one or more input VCFs according to specified intervals
"""


import argparse
import pybedtools as pbt
import pysam
from g2cpy import is_inside, relocate_uri
from os import path


def _open_input_vcf(vpath):
    """
    Open an input VCF, after localizing if necessary
    """

    vbase = path.basename(vpath)

    if path.isfile(vpath):
        return pysam.VariantFile(vpath)

    elif path.isfile(vbase):
        return pysam.VariantFile(vbase)
    
    elif vpath.startswith('gs://'):
        try:
            relocate_uri(vpath, './', verbose=True)
        except:
            raise OSError('Unable to localize {} with gsutil cp'.format(vpath))
        return pysam.VariantFile(vbase)
    
    else:
        raise ValueError('Unable to determine location of {}'.format(vpath))


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vcf-list', metavar='txt', required=True,
                        help='List of input VCFs. Can either be local file paths'
                        ' or GCP URIs, which will be automatically localized. Required.')
    parser.add_argument('--intervals', required=True, metavar='BED',
                        help='BED4 file of desired output shard intervals. ' +
                        'Fourth column must specify output shard VCF prefix. ' +
                        'Required')
    args = parser.parse_args()

    # Read list of VCFs
    with open(args.vcf_list) as fin:
        in_vcfs = list(set([l.rstrip() for l in fin.readlines()]))

    # Determine header for output VCFs based on first input VCF
    header = _open_input_vcf(in_vcfs[0]).header

    # Open connections to output VCFs
    out_map = dict()
    for feature in pbt.BedTool(args.intervals):
        chrom = str(feature.chrom)
        start = int(feature.start)
        end = int(feature.end)
        fname = str(feature[3]) + '.vcf.gz'
        if chrom not in out_map.keys():
            out_map[chrom] = dict()
        out_map[chrom][(start, end, )] = pysam.VariantFile(fname, 'w', header=header)

    # Iterate over each input VCF and route their records accordingly
    for ivcf_path in in_vcfs:
        ivcf = _open_input_vcf(ivcf_path)
        for record in ivcf:
            if record.chrom not in out_map.keys():
                continue
            # Route records based on their position vs. the interval
            for bounds, outvcf in out_map[record.chrom].items():
                if is_inside(record.pos, bounds):
                    outvcf.write(record)
    
    # Close all VCF handles
    for contig in out_map.keys():
        for fh in out_map[contig].values():
            fh.close()


if __name__ == '__main__':
    main()

