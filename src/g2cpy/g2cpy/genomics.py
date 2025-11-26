#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Generic genomics-related helper functions
Note that functions related to VCF parsing or manipulation are found in vcf.py
"""


import pybedtools as pbt
import networkx as nx
import subprocess
from numpy import nanmedian
from re import sub


def bgzip(filename, return_new_fn=False):
    """
    Bgzip a file
    """

    subprocess.run(['bgzip', '-f', filename])

    if return_new_fn:
        return filename + '.gz'


def chrom2int(contig):
    """
    Converts a human primary contig string to an integer for simpler sorting
    """

    contig_order = {'chr' + str(k + 1) : k + 1 for k in range(22)}
    contig_order.update({'chrX' : 23, 'chrY' : 24, 'chrM' : 25})

    return contig_order[contig]


def chromsort(contigs):
    """
    Sorts an iterable of human contigs by their positional order
    """

    return sorted(contigs, key=lambda k: chrom2int(str(k)))


def classify_variant(ref, alt, var_len=None):
    """
    Maps a variant to its G2C variant class & subclass codes
    Returns two values: class, subclass
    """

    pur = 'A G'.split()
    pyr = 'C T'.split()
    nucs = pur + pyr

    # Compute variant length if it is not already defined
    if var_len is None:
        var_len = abs(len(ref) - len(alt))

    # SNVs have zero size and cannot contain multinucleotide refs or alts
    if var_len == 0 and len(ref) == 1 and len(alt) == 1:
        if ref in pur and alt in pur:
            sc = 'ti'
        elif ref in pyr and alt in pyr:
            sc = 'ti'
        else:
            sc = 'tv'
        return 'snv', sc

    # Indels are 1-49bp, excluding certain SV types with indefinite sizes
    elif var_len < 50 and 'CTX' not in alt:
        if 'INS' in alt or 'DUP' in alt:
            return 'indel ins'.split()
        elif 'DEL' in alt:
            return 'indel del'.split()
        elif len(ref) > len(alt):
            return 'indel del'.split()
        elif len(ref) < len(alt):
            return 'indel ins'.split()
        elif len(ref) == len(alt):
            return 'indel cpx'.split()
        else:
            msg = 'g2cpy classify_variant unable to classify suspected indel ' + \
                  'with ref {}, alt {}, and length {:,}'
            exit(msg.format(ref, alt, var_len))

    # SVs are >=50bp, except for rare types with indefinite sizes
    else:
        if alt.startswith('<') and alt.endswith('>'):
            sc = sub('<|>', '', alt).split(':')[0]
            return 'sv', sc
        elif len(ref) > len(alt):
            return 'sv DEL'.split()
        elif len(ref) < len(alt):
            return 'sv INS'.split()
        else:
            msg = 'g2cpy classify_variant unable to classify suspected SV ' + \
                  'with ref {}, alt {}, and length {:,}'
            exit(msg.format(ref, alt, var_len))


def cluster_bedtool(bt, recip=10e-10):
    """
    Performs single-linkage clustering and merging of a pbt.BedTool
    Takes an `r` parameter to enforce reciprocal overlap between elements before merging
    """

    # Format intervals to be clustered so each row has a unique identifier
    str_feats = '\n'.join(['{}\t{}\t{}\trow_{}'.format(*f, i) 
                           for i, f in enumerate(bt.cut(range(3)).sort())])
    bt = pbt.BedTool(str_feats, from_string=True)

    # Compute hits per `recip`
    hits = bt.intersect(bt, r=True, f=recip, wa=True, wb=True)

    # Build graph of qualifying hits
    g = nx.Graph()
    for hit in hits:
        g.add_node(hit[3], chrom=hit[0], start=int(hit[1]), end=int(hit[2]))
        if hit[3] != hit[7]:
            g.add_node(hit[7], chrom=hit[4], start=int(hit[5]), end=int(hit[6]))
            g.add_edge(hit[3], hit[7])

    # Iterate over subgraphs and compute median coordinates for each
    out_feats = []
    for cluster in nx.connected_components(g):
        chrom = g.nodes[list(cluster)[0]]['chrom']
        pos = int(nanmedian([g.nodes[n]['start'] for n in cluster]))
        end = int(nanmedian([g.nodes[n]['end'] for n in cluster]))
        out_feats.append([chrom, str(pos), str(end)])
    out_feat_s = ['\t'.join(v) for v in 
                  sorted(out_feats, key=lambda x: (x[1], x[2]))]
    
    return pbt.BedTool('\n'.join(out_feat_s), from_string=True)


def determine_filetype(path, return_extension=False):
    """
    Determine file extension for common genomic data formats
    """

    # Enumerate candidate suffix matches
    suf_dict = {'cram' : ['cram'],
                'bam' : ['bam'],
                'vcf' : ['vcf'],
                'compressed-vcf' : 'vcf.gz vcf.bgz vcf.gzip vcf.bgzip'.split(),
                'bed' : ['bed'],
                'compressed-bed' : 'bed.gz bed.bgz bed.gzip bed.bgzip'.split(),
                'bigwig' : '.bw .bigwig .bigWig .BigWig'.split(),
                'fasta' : '.fa .fasta'.split(),
                'compressed-fasta' : '.fa.gz .fa.gzip .fasta.gz .fasta.gzip'.split(),
                'hic' : ['hic'],
                'gtf' : ['gtf'],
                'compressed-gtf' : ['gtf.gz'],
                'gatk-interval-list' : '.intervals',
                'picard-interval-list' : '.interval_list',
                'tsv' : ['.tsv'],
                'compressed-tsv' : ['.tsv.gz']}

    for ftype, suffs in suf_dict.items():
        suf_hits = [s for s in suffs if path.endswith(s)]
        if len(suf_hits) > 0:
            if return_extension:
                return ftype, suf_hits[0]
            else:
                return ftype

    # If no matches are found, return None
    if return_extension:
        return None, None
    else:
        return None


def make_tabix_index(filename):
    """
    Index a bgzipped file with tabix
    """

    subprocess.run(['tabix', '-f', filename])


def name_variant(chrom, pos, ref, alt, vc, vsc, varlen, suffix_length=10):
    """
    Determine unique G2C variant identifier based on its site information
    """

    if vc == 'snv' \
    or (vc == 'indel' and (len(ref) + len(alt) + 1) <= suffix_length):
        return '_'.join([str(x) for x in [chrom, pos, ref, alt]])
    else:
        vhash_str = ''.join([str(x) for x in [vsc, pos, ref, varlen, alt]])
        vhash = hash_string(vhash_str, out_length=suffix_length)
        return '_'.join([str(x) for x in [chrom, pos, vsc, vhash]])

