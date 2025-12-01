#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Divide a GATK-style interval_list file into shards, either yielding shards matched
on a fixed target genomic size or balanced for the density of variants in a reference
dataset, like gnomAD
"""


# Import libraries
import argparse
import gzip
import numpy as np
import pandas as pd
import pybedtools as pbt
import re
from g2cpy import chrom2int, determine_filetype
from sys import stdin, stdout
from pysam import TabixFile


def uniform_split(int_list, n_shards, min_size=5000, zero_based=False, verbose=False):
    """
    Split each interval into the same number of equal-sized shards such that the
    total number of shards matches n_shards
    """

    # Assign buffers for size and start based on coordinate system
    if zero_based:
        coord_buf = 0
    else:
        coord_buf = 1

    # Sort intervals by size, from largest to smallest
    int_list = sorted(int_list, key=lambda x: x[0], reverse=True)

    # Determine sharding parameters
    target_spi = int(np.floor(n_shards / len(int_list)))
    max_spi = [int(np.floor(size / min_size)) for size, vals in int_list]
    spi = [np.nanmin([target_spi, m]) for m in max_spi]
    rem = n_shards - sum(spi)
    while rem > 0:
        for i in range(len(spi)):
            if spi[i] < max_spi[i]:
                spi[i] += 1
                rem -= 1
                if rem == 0:
                    break

    # Uniformly shard each interval
    shards = []
    for i, ivals in enumerate(int_list):
        isize = ivals[0]
        ifields = ivals[1]
        ishards = spi[i]

        chrom = ifields[0]
        start = int(ifields[1])
        end = int(ifields[2])

        if ishards < 2:
            shards.append([chrom, start, end])
            continue

        split_points = np.round(np.linspace(start, end, ishards + 1), 0)
        for si in range(ishards):
            cstart = int(split_points[si])
            if si > 0:
                cstart += coord_buf
            cend = int(split_points[si+1])
            csize = cend - cstart + coord_buf
            shards.append([chrom, cstart, cend])

    return sorted(shards, key=lambda x: (chrom2int(x[0]), int(x[1]), int(x[2])))


def split_by_density(int_list, var_beds, vars_per_shard, min_size=5000, 
                     max_size=np.inf, zero_based=False, verbose=False):
    """
    Splits intervals in int_list into shards where each shard carries no more than
    vars_per_shard total variants present in var_beds
    """

    # Convert intervals to pbt.BedTool and ensure coordinate-sorted
    int_bt = pbt.BedTool('\n'.join(int_list), from_string = True).sort()

    # Open tabix handles to var_beds
    var_tbx = {i : TabixFile(var_beds[i]) for i in range(len(var_beds))}

    # Iterates over each original input interval
    shards = []
    for orig_int in int_bt:
        chrom = orig_int.chrom
        istart = orig_int.start
        iend = orig_int.end

        # Open handles to var beds for positions within this original interval
        ivar_tbxs = {i : var_tbx[i].fetch(chrom, istart, iend) for i, b in var_tbx.items()}

        # Initialize coordinate pointers for this shard
        def _increment_tbx(tbx):
            try:
                line = next(tbx)
                return int(line.split('\t')[1])
                # Once the tabix handle is exhausted, emit an impossibly large coordinate
            except:
                return 10e10
        tbx_pointers = {i : _increment_tbx(t) for i, t in ivar_tbxs.items()}

        # Function to consume current minimum position, increment that pointer, and return updated tbx_pointers
        def _spend_variant(tbx_pointers):
            i = np.argmin(list(tbx_pointers.values()))
            tbx_pointers[i] = _increment_tbx(ivar_tbxs[i])
            return tbx_pointers

        # Count variants until either vars_per_shard have been allocated
        # or until iend (or max_size) is reached, in which case always emit the
        # shard and start a new one
        s_start = istart
        cur_pos = np.min(list(tbx_pointers.values()) + [iend])
        k = 0
        s_size = 0
        while cur_pos <= iend:

            # If current position is sufficiently close to end of interval s/t
            # a new shard larger than min_size cannot fit between cur_pos and iend,
            # round cur_pos up to iend and break this loop
            if cur_pos >= iend - min_size:
                cur_pos = iend
                break

            # Find next variant, increment that pointer, add to tally, and update shard size
            tbx_pointers = _spend_variant(tbx_pointers)
            k += 1
            cur_pos = np.min(list(tbx_pointers.values()) + [iend])
            s_size = cur_pos - s_start
            if not zero_based:
                s_size += 1

            # Once the number of variants has been tallied or max_size has been 
            # reached, emit current shard and update new shard start/end
            if k >= vars_per_shard or s_size >= max_size:
                shards.append([chrom, s_start, cur_pos])
                if verbose:
                    msg = 'Split interval {}:{:,}-{:,} at {}:{:,}-{:,} ' + \
                          '({:.0f} kb) after {:,} variants\n'
                    stdout.write(msg.format(chrom, istart, iend, chrom, s_start, 
                                            cur_pos, (cur_pos - s_start) / 1000, k))
                if not zero_based:
                    cur_pos += 1
                s_start = cur_pos
                cur_pos += 1
                k = 0
                s_size = 0

        # As the above loop should always break with imperfect sharding, we
        # always need to emit the final shard to avoid right-truncation
        shards.append([chrom, s_start, int(np.max([cur_pos, iend]))])
        if verbose:
            msg = 'Yielding {}:{:,}-{:,} ({:.0f} kb) as the remainder of ' + \
                  'interval {}:{:,}-{:,} after {:,} variants\n'
            stdout.write(msg.format(chrom, s_start, cur_pos, (cur_pos - s_start) / 1000, 
                                    chrom, istart, iend, k))


    return sorted(shards, key=lambda x: (chrom2int(x[0]), int(x[1]), int(x[2])))


def split_by_size(int_list, target_size, max_size=np.inf, zero_based=False):
    """
    Recursively splits the largest interval in half until mean interval is ≤ --target-size
    Also ensures no intervals are larger than max_size; such intervals will be 
    split recursively until all shards are < max_size
    """

    # Assign buffers for size and start based on coordinate system
    if zero_based:
        coord_buf = 0
    else:
        coord_buf = 1

    int_df = pd.DataFrame(int_list, columns='isize fields'.split())

    # Take a first pass through int_df and divide all intervals > max_size into 
    # their component parts
    split_idx = np.where(int_df.isize > max_size)[0]
    children = []
    for big_idx in split_idx:

        # Get info for parent interval
        parent = int_df.fields[big_idx]
        chrom = parent[0]
        parent[1] = int(parent[1])
        parent[2] = int(parent[2])

        # Split parent into the smallest number of uniform children 
        # s/t no child is > max_size
        n_splits = int(np.ceil((parent[2] - parent[1] + coord_buf) / max_size))
        split_points = np.round(np.linspace(parent[1], parent[2], n_splits + 1), 0)
        for si in range(n_splits):
            cstart = int(split_points[si])
            if si > 0:
                cstart += coord_buf
            cend = int(split_points[si+1])
            csize = cend - cstart + coord_buf
            children.append((csize, [chrom, cstart, cend, *parent[3:]]))

    # Drop supersized intervals and append their split children
    if len(split_idx) > 0:
        int_df = pd.concat([int_df.drop(split_idx),
                    pd.DataFrame(children, columns=int_df.columns)]).\
            reset_index(drop=True)

    # After splitting all supersized intervals, continue recursively splitting 
    # the largest remaining interval in half until the average shard size is 
    # <= the target size
    avg_size = int_df.isize.mean()
    while avg_size > target_size:

        # Find largest interval
        big_idx = int_df.isize.argmax()
        parent = int_df.fields[big_idx]
        parent[1] = int(parent[1])
        parent[2] = int(parent[2])

        # Divide interval in half
        midpoint = int(np.floor(np.mean([parent[2], parent[1]])))
        children = [(midpoint - parent[1], [*parent[0:2], midpoint, *parent[3:]], ),
                    (parent[2] - midpoint, [parent[0], midpoint+1, *parent[2:]], ),]

        # Drop parent and append children
        int_df = pd.concat([int_df.drop(big_idx),
                            pd.DataFrame(children, columns=int_df.columns)]).\
                    reset_index(drop=True)

        # Recompute mean interval size
        avg_size = int_df.isize.mean()

    # Return list of new intervals sorted by start coordinate
    return sorted(int_df.fields.tolist(), 
                  key=lambda x: (chrom2int(x[0]), int(x[1]), int(x[2])))


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input-intervals', required=True, 
                        help='Path to input .interval_list', default='stdin')
    parser.add_argument('-t', '--target-size', help='Desired average interval ' +
                        'size. By default, intervals will be divided based on ' +
                        'their genomic size, controlled by this parameter.',
                        type=float, default=10e10)
    parser.add_argument('--n-shards', type=int, help='Total number of desired ' +
                        'shards. Specifying this option will override ' +
                        '--target-size, but the exact behavior depends on other ' +
                        'options. If --vars-per-shard and --var-sites are also ' +
                        'specified, will perform density-based splitting. ' +
                        'Otherwise, will uniformly divide each interval the same ' +
                        'number of times such that --n-shards in total are produced.')
    parser.add_argument('--var-sites', action='append', help='One or more BED ' +
                        'files listing known variant sites in an external ' +
                        'reference dataset, like gnomAD. May be provided ' +
                        'multiple times. If specified with --vars-per-shard, ' +
                        'intervals will be divided into shards balanced by ' +
                        'number of variants per shard.')
    parser.add_argument('--vars-per-shard', type=float, help='Desired number ' +
                        'of variants present in --var-sites to allocate per shard. ' +
                        'Required for variant density-based shard balancing ' +
                        'if --var-sites is also optioned.')
    parser.add_argument('-d', '--maximum-distance-between-intervals', type=int,
                        help='Maximum distance to permit between adjacent ' +
                        'intervals before writing to a new output file. By ' +
                        'default, all output intervals will be written to one file')
    parser.add_argument('--min-interval-size', type=int, default=5000,
                        help='Size of smallest interval shard to be emitted, in ' +
                        'base pairs.')
    parser.add_argument('--max-interval-size', type=int, default=np.inf,
                        help='Size of largest interval shard to be emitted, ' +
                        'in base pairs. [default: no maximum]')
    parser.add_argument('--gatk-style', action='store_true', default=False,
                        help='Write output as simple GATK-style intervals ' +
                        '[default: write Picard-style intervals with header]')
    parser.add_argument('--bed-style', action='store_true', default=False,
                        help='Write output as BED intervals with 0-based coordinates ' +
                        '[default: write Picard-style intervals with header]')
    parser.add_argument('-o', '--output-intervals', type=str, default='stdout',
                        help='Path to output .interval_list')
    parser.add_argument('-p', '--output-prefix', type=str,
                        help='Filename prefix for interval_list shards. Only ' +
                        'used if --maximum-distance-between-intervals is ' +
                        'specified, otherwise all intervals will be written to ' +
                        '--output-intervals')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print ' +
                        'diagnostic logging to stdout')
    args = parser.parse_args()

    # Determine & report desired splitting mode 
    if args.n_shards:
        if args.vars_per_shard is not None \
        and len(args.var_sites) > 0:
            split_mode = 'density'
            msg = 'Splitting intervals by variant density due to the presence of ' + \
                  '--n-shards, --var-sites, and --vars-per-shard'
        else:
            split_mode = 'uniform'
            msg = 'Uniformly splitting each interval the same number of times due ' + \
                  'to the presence of --n-shards on its own'
    else:
        split_mode = 'size'
        msg = 'Splitting intervals by size due to the lack of other specific options'
    print(msg + '\n')

    # Open connection to input file
    if args.input_intervals in '- stdin /dev/stdin'.split():
        infile = stdin
    else:
        if 'compressed' in determine_filetype(args.input_intervals):
            infile = gzip.open(args.input_intervals, 'rt')
        else:
            infile = open(args.input_intervals)

    # Open connection to output file
    multi_out = args.maximum_distance_between_intervals is not None
    if multi_out:
        if not args.gatk_style:
            outfile_suffix = 'interval_list'
        else:
            outfile_suffix = 'intervals'
        outfile_fmt = args.output_prefix + '.{}.' + outfile_suffix
        n_outfiles = 0
        outfile = None
    else:
        if args.output_intervals in '- stdout /dev/stdout'.split():
            outfile = stdout
        else:
            outfile = open(args.output_intervals, 'w')

    # Read input file into memory
    int_list = []
    if multi_out:
            out_header = ''
    for line in infile.readlines():

        # Write header from input to output
        if line.startswith('@'):
            if not args.gatk_style:
                if multi_out:
                    out_header += line
                else:
                    outfile.write(line)

        # Store input intervals as tuples of (interval size, list of fields)
        # or as strings for pbt.BedTool, depending on split_mode
        else:
            if '\t' in line:
                fields = line.rstrip().split()
            else:
                fields = re.split('[:-]', line.rstrip())
            if split_mode == 'density':
                int_list.append('\t'.join(fields[:3]))
            else:
                size = int(fields[2]) - int(fields[1])
                int_list.append((size, fields, ))

    # After reading input intervals, perform splitting based on value of split_mode
    if split_mode == 'uniform':
        shards = uniform_split(int_list, args.n_shards, args.min_interval_size,
                               args.bed_style, args.verbose)
    elif split_mode == 'density':
        shards = split_by_density(int_list, args.var_sites, args.vars_per_shard, 
                                  args.min_interval_size, args.max_interval_size,
                                  args.bed_style, args.verbose)
    else:
        shards = split_by_size(int_list, args.target_size, args.max_interval_size,
                               args.bed_style)

    # Sort sharded intervals and write to output file
    prev_chrom = 'notreal'
    prev_end = -9e8
    for vals in sorted(shards, key=lambda x: (chrom2int(x[0]), int(x[1]), int(x[2]))):
        # Format output line
        if args.gatk_style:
            outline = '{}:{}-{}'.format(*vals[:3])
        elif args.bed_style:
            outline = '{}\t{}\t{}'.format(*vals[:3])
        else:
            outline = '\t'.join([str(k) for k in vals[:3] + ['+', '. intersection ACGTmer']])

        # Write to output file depending on desired output behavior
        if multi_out:
            # Check if a new outfile is necessary
            if vals[0] != prev_chrom \
            or int(vals[1]) > prev_end + args.maximum_distance_between_intervals:
                if outfile is not None:
                    outfile.close()
                n_outfiles += 1
                outfile = open(outfile_fmt.format(n_outfiles), 'w')
                outfile.write(out_header)
        outfile.write(outline + '\n')
        prev_chrom = vals[0]
        prev_end = int(vals[2])

    # Clear buffer
    outfile.close()


if __name__ == '__main__':
    main()

