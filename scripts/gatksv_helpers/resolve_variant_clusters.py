#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Extract and collapse clusters of records from an input VCF
"""


import argparse
import csv
import g2cpy
import numpy as np
import pysam
from sys import stdin, stdout


def write_nonredundant_records(records, outvcf):
    """
    Given a list of pysam.VariantRecord objects, write those with 
    unique POS, END, and ALT values to outvcf
    """

    pe = set([(r.pos, r.stop) for r in records])
    for pos, end in sorted(list(pe)):
        for rec in records:
            if rec.pos == pos and rec.stop == end:
                outvcf.write(rec)
                break


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--in-vcf', metavar='vcf', default='stdin',
                        help='input .vcf [default: stdin]')
    parser.add_argument('--clusters', metavar='bed', help='BED4 with one line ' +
                        'per cluster to process. Fourth column is comma-delimited ' +
                        'list of cluster members', required=True)
    parser.add_argument('--out-vcf', metavar='vcf', default='stdout',
                        help='output .vcf [default: stdout]')
    parser.add_argument('-p', '--prefix', help='Name prefix for clustered variants',
                        default='reclustered')
    parser.add_argument('-b', '--buffer', type=int, default=10, 
                        help='Buffer for VCF queries, in base pairs')
    args = parser.parse_args()

    # Open connection to input vcf
    if args.in_vcf in 'stdin /dev/stdin -'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.in_vcf)

    # Open connection to cluster BED file
    cluster_bed = csv.reader(open(args.clusters), delimiter='\t')

    # Get list of samples from input vcf
    samples = [s for s in invcf.header.samples]

    # Open connection to output vcf
    if args.out_vcf in 'stdout /dev/stdout -'.split():
        outvcf = pysam.VariantFile(stdout, mode='w', header=invcf.header)
    else:
        outvcf = pysam.VariantFile(args.out_vcf, mode='w', header=invcf.header)

    # Process each cluster in serial
    k = 0
    for chrom, start, end, mem_str in cluster_bed:
        
        vids = mem_str.split(',')

        # Gather records from VCF
        records = []
        for rec in invcf.fetch(chrom, 
                               np.nanmax([int(start) - args.buffer, 0]), 
                               int(end) + args.buffer):
            if rec.id in vids:
                records.append(rec)

        # Check that all records can be found
        if len(records) != len(vids):
            msg = 'Unable to find all {:,} records in --in-vcf for cluster at {}:{}-{}'
            exit(msg.format(len(vids), chrom, start, end))
        k += 1

        # Get consensus coordinates of merged record
        cpos = int(np.floor(np.nanmedian([r.pos for r in records])))
        cend = int(np.floor(np.nanmedian([r.stop for r in records])))

        # For complex variants, further ensure all variants have the same
        # subtype (INFO/CPX_TYPE) and their CPX_INTERVALS can be properly 
        # unified before merging into a single record
        if records[0].info.get('SVTYPE') == 'CPX':
            cpx_types = list(set([r.info.get('CPX_TYPE') for r in records]))
            if len(cpx_types) > 1:
                write_nonredundant_records(records, outvcf)
                continue
            try:
                cpx_ints = g2cpy.integrate_cpx_intervals(records, cpx_types[0], 
                                                         chrom, cpos, cend)
            except:
                write_nonredundant_records(records, outvcf)
                continue

        # Use first record as a template
        newrec = records[0].copy()

        # Assign basic (non-INFO) record information
        newrec.pos = cpos
        newrec.stop = cend
        newrec.id = '{}_{}'.format(args.prefix, k)
        newrec.qual = int(np.round(np.nanmean([r.qual for r in records])))
        newrec.filter.clear()
        for f in list(set(g2cpy.recursive_flatten([r.filter.items() for r in records]))):
            newrec.filter.add(f)

        # Merge INFOs
        newrec.info.clear()
        newrec.info.update(g2cpy.integrate_infos(records, do_cpx_intervals=False))
        if newrec.info['SVTYPE'] == 'CPX':
            newrec.info['CPX_INTERVALS'] = cpx_ints
        if newrec.info.get('SVLEN', 0) < 50:
            newrec.info['SVLEN'] = int(np.nanmax([newrec.stop - newrec.pos, 0]))
        if newrec.info['SVLEN'] == 'INS':
            import pdb; pdb.set_trace()

        # Merge GTs
        newrec = g2cpy.integrate_gts(newrec, records)

        # Update AC/AN/AF
        af_stats = g2cpy.compute_allele_freq_stats(newrec, keys='AN AC AF'.split())
        for key, value in af_stats.items():
            newrec.info[key] = value

        # Write clustered record to output VCF
        outvcf.write(newrec)

    # Clear buffers
    outvcf.close()


if __name__ == '__main__':
    main()

