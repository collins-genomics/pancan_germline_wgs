#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Identify and tag or correct variants in a multi-sample VCF with markedly different
frequencies between user-specified subsets of samples (e.g., batches, cohorts, etc.)
"""


import argparse
import numpy as np
import pysam


def build_group_map(group_membership_tsv, all_samples):
    """
    Load a group membership map from .tsv as dict
    """

    group_map = {}
    group_members = {}

    with open(group_membership_tsv) as map_in:
        for line in map_in.readlines():
            gid, sid = line.rstrip().split('\t')

            # Only add samples present in all_samples (taken from VCF header)
            if sid not in all_samples:
                continue

            if sid not in group_map.keys():
                group_map[sid] = set()
            group_map[sid].add(gid)
            if gid not in group_members:
                group_members[gid] = set()
            group_members[gid].add(sid)

    return group_map, group_members


def is_multiallelic(record):
    """
    Check if record is multiallelic
    """

    if len(record.alts) > 1:
        return True

    mcnv_alts = '<CNV> <MCNV>'.split()
    if any([a in mcnv_alts for a in record.alts]):
        return True

    mcnv_svtypes = 'CNV MCNV'.split()
    if record.info.get('SVTYPE') in mcnv_svtypes:
        return True

    if 'MULTIALLELIC' in record.filter.keys() \
    or 'MULTIALLELIC' in record.info.keys():
        return True

    return False


def collect_batch_effect_info(record, group_map):
    """
    Collect information required for batch effect detection
    """

    bfx_info = {}

    for sid, sinfo in record.samples.items():
        
        # Only count non-null GTs
        a = [a for a in sinfo.get('GT', (None, None)) if a is not None]
        if len(a) == 0:
            continue

        gids = group_map.get(sid)
        
        if gids is None:
            continue

        if len(gids) == 0:
            continue

        for gid in gids:
            if gid not in bfx_info.keys():
                bfx_info[gid] = {'total' : 0, 'nonref' : 0}
            bfx_info[gid]['total'] += 1
            if any([k > 0 for k in a]):
                bfx_info[gid]['nonref'] += 1

    for gid in bfx_info.keys():
        if bfx_info[gid]['total'] == 0:
            bfx_info[gid]['frac'] = np.nan
        else:
            bfx_info[gid]['frac'] = bfx_info[gid]['nonref'] / bfx_info[gid]['total']

    return bfx_info


def check_batch_effects(bfx_info, min_n, lower_freq, upper_freq):
    """
    Identify which groups meet criteria as either unusually low or high frequency
    """

    low_groups = [gid for gid, gdat in bfx_info.items() \
                  if gdat['total'] >= min_n \
                  and gdat['frac'] < lower_freq]

    high_groups = [gid for gid, gdat in bfx_info.items() \
                   if gdat['total'] >= min_n \
                   and gdat['frac'] > upper_freq]

    any_bfx = len(low_groups) > 0 and len(high_groups) > 0

    return any_bfx, low_groups, high_groups


def pool_freqs(bfx_info, groups):
    """
    Pool frequencies for one or more groups
    """

    total, nonref = 0, 0
    for gid in groups:
        total += bfx_info[gid]['total']
        nonref += bfx_info[gid]['nonref']
    
    return {'total' : total, 'nonref' : nonref, 'frac' : nonref / total}


def mask_gts(record, group_members, groups_to_mask, iflag):
    """
    Mask all genotypes for samples from one or more specified groups_to_mask
    """

    for gid in groups_to_mask:
        for sid in group_members[gid]:
            gt = record.samples[sid]['GT']
            record.samples[sid]['GT'] = tuple([None] * len(gt))

    record.info[iflag] = True
    
    return record


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input-vcf', required=True, metavar='VCF', 
                        type=str, help='Input .vcf. Required.')
    parser.add_argument('-m', '--group-membership', required=True, type=str,
                        help='Two-column .tsv mapping group name (column 1) ' +
                        'to sample ID (column 2). One row per sample. Required.',
                        metavar='tsv')
    parser.add_argument('-o', '--output-vcf', help='Output .vcf [default: stdout]',
                        default='stdout', metavar='VCF', type=str)
    parser.add_argument('-S', '--strict', default=False, action='store_true',
                        help='Operate in strict mode, which tags records ' +
                        'with a non-PASS FILTER if any qualifying ' +
                        'inter-group differences are observed')
    parser.add_argument('--filter-id', type=str, metavar='String',
                        help='Custom VCF FILTER ID to use when tagging records ' +
                        'with inter-group differences. Only used with --strict.')
    parser.add_argument('--filter-description', type=str, metavar='String',
                        help='Custom VCF FILTER description to use when ' +
                        'tagging records with inter-group differences. Only ' +
                        'used with --strict.')
    parser.add_argument('--info-flag-id', type=str, metavar='String',
                        help='Custom VCF INFO ID to use when flagging records ' +
                        'with GTs masked due to inter-group differences.',
                        default='GROUPWISE_MASKED_GENOTYPES')
    parser.add_argument('-M', '--min-samples', default=10, type=int, 
                        help='Minimum number of samples per group to be ' +
                        'included in comparisons [default: 10]', metavar='Int')
    parser.add_argument('-L', '--lower-freq', default=0.1, type=float,
                        help='Lower frequency threshold for inter-group ' +
                        'comparisons [default: 0.1]', metavar='Float')
    parser.add_argument('-U', '--upper-freq', default=0.5, type=float,
                        help='Upper frequency threshold for inter-group ' +
                        'comparisons [default: 0.5]', metavar='Float')
    args = parser.parse_args()

    # Open connection to input VCF
    invcf = pysam.VariantFile(args.input_vcf)
    all_samples = set([s for s in invcf.header.samples])

    # Read map of group memberships as dict
    group_map, group_members = build_group_map(args.group_membership, all_samples)

    # Ensure at least --min-samples from at least two groups are present in the input VCF header
    if len(group_members.keys()) < 2:
        errmsg = 'At least two different groups must be present in ' + \
                 'the first column of --group-membership. Exiting.'
        exit(errmsg)
    n_per_group = [len(gsamps.intersection(all_samples)) for gsamps in group_members.values()]
    if len([k for k in n_per_group if k >= args.min_samples]) < 2:
        errmsg = 'At least --min-samples ({}) samples from at least two ' + \
                 'different groups in --group-membership must be present ' + \
                 'in --input-vcf. Exiting.'
        exit(errmsg)

    # Modify output VCF header to have necessary tags depending on mode
    out_header = invcf.header
    i_line = '##INFO=<ID={},Number=0,Type=Flag,' + \
             'Description="A subset of GTs for this record have been ' + \
             'masked by post hoc application of mark_batch_effects.py">'
    if args.info_flag_id not in out_header.info.keys():
        out_header.add_line(i_line.format(args.info_flag_id))
    if args.strict:
        filter_id = 'NONUNIFORM_GENOTYPES'
        filter_descrip = 'Variant exhibits excessive differences in ' + \
                         'non-ref GT rates between prespecified subsets ' + \
                         'of samples (e.g., batches)'
        if args.filter_id is not None:
            filter_id = args.filter_id
        if args.filter_description is not None:
            filter_descrip = args.filter_description
        f_line = '##FILTER=<ID={},Description="{}">'
        if filter_id not in out_header.filters.keys():
            out_header.add_line(f_line)
        out_header.add_line(f_line.format(filter_id, filter_descrip))

    # Open connection to output VCF
    if args.output_vcf in 'stdout /dev/stdout -'.split():
        from sys import stdout
        outvcf = pysam.VariantFile(stdout, 'w', header=out_header)
    else:
        outvcf = pysam.VariantFile(args.output_vcf, 'w', header=out_header)

    # Process each sample in serial and write to --output-vcf
    for record in invcf.fetch():

        # Skip multiallelic records
        if is_multiallelic(record):
            outvcf.write(record)
            continue

        bfx_info = collect_batch_effect_info(record, group_map)
        any_bfx, low_groups, high_groups = \
            check_batch_effects(bfx_info, args.min_samples, 
                                args.lower_freq, args.upper_freq)

        # Skip records with no observed groupwise differences
        if not any_bfx:
            outvcf.write(record)
            continue

        # Mask GTs from batches with the most aberrant non-ref GT rates
        # Note that the "correct" set of groups is determined based on average
        # distance to the full sample frequency
        all_frac = pool_freqs(bfx_info, bfx_info.keys())['frac']
        low_frac = pool_freqs(bfx_info, low_groups)['frac']
        high_frac = pool_freqs(bfx_info, high_groups)['frac']
        low_dist = abs(all_frac - low_frac)
        high_dist = abs(all_frac - high_frac)
        if low_dist > high_dist:
            bad_groups = low_groups
        else:
            bad_groups = high_groups
        record = mask_gts(record, group_members, bad_groups, args.info_flag_id)

        # Add non-PASS FILTER if executed in --strict mode
        if args.strict:
            record.filter.add(filter_id)

        outvcf.write(record)

    # Close connection to output VCF to clear buffer
    outvcf.close()


if __name__ == '__main__':
    main()

