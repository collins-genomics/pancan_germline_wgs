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


common_numeric_feats = 'GQ DP MIN_DP'.split()
gatksv_numeric_feats = 'CNQ OGQ PE_GT PE_GQ RD_GQ SL SR_GT SR_GQ'.split()



def gather_gt_features(sdat, vc, need_header=False):
    """
    Main function to collect genotype-level features
    """

    outvals, outcols = [], []

    # Collect common numeric features defined for all VCs
    for nf in common_numeric_feats:
        outvals.append(sdat.get(nf, '.'))
    if need_header:
        outcols += common_numeric_feats

    # Collect ref PL from GATK-HC as a ratio vs. max(PL)
    ref_pl_ratio = '.'
    if 'PL' in sdat.keys():
        pls = list(sdat.get('PL', (None, None, None, )))
        if all([x is not None for x in pls]):
            ref_pl = pls[0]
            max_pl = max(pls)
            ref_pl_ratio = ref_pl / max_pl
    outvals.append(ref_pl_ratio)
    if need_header:
        outcols.append('ref_pl_ratio')

    # Collect numeric features specific to GATK-SV
    if vc in 'indel sv'.split():
        for nf in gatksv_numeric_feats:
            outvals.append(sdat.get(nf, '.'))
        if need_header:
            outcols += gatksv_numeric_feats

        # Absolute copy number and marginal difference vs. expected ploidy
        cn = sdat.get('RD_CN', None)
        outvals.append(cn)
        ecn = sdat.get('ECN', None)
        dcn = None
        if cn is not None and ecn is not None:
            try:
                dcn = int(cn) - int(ecn)
            except:
                pass
        outvals.append(dcn)
        if need_header:
            outcols += 'rd_cn dcn'.split()

        # Store EV as sorted string for downstream tokenizing
        evl = [x for x in sdat.get('EV', tuple()) if x is not None]
        if len(evl) == 0:
            ev = '.'
        else:
            ev = ','.join(sorted(list(evl)))
        outvals.append(ev)
        if need_header:
            outcols.append('gt_ev')

    # Convert all Nones in outvals
    for i, v in enumerate(outvals):
        if v is None:
            outvals[i] = '.'

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
                        type=str, help='Output .tsv. First two columns are ' +
                        'chrom and pos, so can be block compressed and indexed')
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
    if args.output_tsv in 'stdout /dev/stdout -'.split():
        outfile = stdout
    else:
        if 'compressed' in g2cpy.determine_filetype(args.output_tsv):
            outfile = gzip.open(args.output_tsv, 'wt', encoding='utf-8')
        else:
            outfile = open(args.output_tsv, 'w')

    # Iterate over input VCF records and samples within each record
    header_written = False
    for i, record in enumerate(invcf):
        vc, vsc = g2cpy.classify_record(record)
        for sid, sdat in record.samples.items():

            # Check if this sample needs to be processed
            sgt = g2cpy.parse_gt(sdat['GT'])
            if sgt['AC'] == 0 and not args.collect_all:
                continue

            # Compute allele balance
            # We prefer to use GATK AD values where available
            # Otherwise, we fall back on the ratio of genotyped alleles
            sac = sgt['AC']
            if 'AD' in sdat.keys() \
            and any([x is not None for x in sdat.get('AD', (None,))]):
                total_ad = int(sum(sdat['AD']))
                alt_ad = int(sum(sdat['AD'][1:]))
                if total_ad > 0:
                    ab = alt_ad / total_ad
            else:
                san = sgt['AN'] - sgt['N_missing']
                if san > 0:
                    ab = sac / san
                else:
                    ab = '.'
            
            # Collect basic info
            outvals = [record.chrom, record.pos, record.id, sid, sac, ab]
            if not header_written:
                outcols = '#chrom pos vid sid sac ab'.split()

            # Collect other features
            add_outvals, add_outcols = gather_gt_features(sdat, vc, not header_written)
            outvals += add_outvals
            outcols += add_outcols

            # Write header if first line
            if not header_written:
                outfile.write('\t'.join([x.lower() for x in outcols]) + '\n')
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

