#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Genomics-related helper functions
"""


import subprocess
import numpy as np
from collections.abc import Iterable
from re import sub
from statistics import multimode
from .gatksv import integrate_cpx_intervals
from .utilities import hash_string, recursive_flatten


def apply_across_samples(record, function=None, field='GT', samples=None):
    """
    Apply a user-specified function across one FORMAT field for all samples in pysam.VariantRecord
    If 'fuction' is None (default), returns all values as a list
    Provide a list of sample IDs as 'samples' to restrict operation to a subset of samples in VCF
    """

    if samples is None:
        samples = list(record.samples.keys())

    vals = [record.samples[sid][field] for sid in samples]

    if function is None:
        return vals
    else:
        return list(map(function, vals))


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


def compute_allele_freq_stats(record, keys=['AC', 'AN', 'AF']):
    """
    Helper wrapper around apply_across_samples() to compute common frequency 
    stats for a pysam.VariantRecord object
    Returns: dict of floats keyed by 'keys'
    """

    res = dict()

    if 'AC' in keys or 'AF' in keys:
        ac = np.nansum(apply_across_samples(record, 
                       lambda x: np.nansum([a > 0 for a in x if a is not None])))
        res['AC'] = ac
    if 'AN' in keys or 'AF' in keys:
        an = np.nansum(apply_across_samples(record, 
                       lambda x: np.nansum([a is not None for a in x])))
        res['AN'] = an
    if 'AF' in keys:
        af = ac / an
        res['AF'] = af

    for ak in res.keys():
        if ak not in keys:
            res.pop(ak)

    return res


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


def integrate_gts(target_record, records, sort_key='GQ'):
    """
    Merge genotypes across pysam.VariantRecords
    GT selection is prioritized as:
    1. Null / no-call (./.)
    2. All non-ref GTs
    3. All ref GTs
    Ties are broken by sort_key, and entire FORMAT is retained from winner GT
    Writes new genotypes into target_record
    """

    if len(records) < 2:
        return target_record

    # Check to ensure all samples are the same across all records
    sids_per_rec = [set(r.samples.keys()) for r in records]
    all_sids = set(recursive_flatten(sids_per_rec))
    if len(set(map(len, [all_sids] + sids_per_rec))) > 1:
        msg = 'Attempted to integrate genotypes but failed due to inconsistent ' + \
              'samples for records {}'
        exit(msg.format(', '.join([r.id for r in records])))

    # Direct access to target GT fields
    newgts = target_record.samples

    def __sort_gts(gt, tiebreak='GQ'):
        """
        Custom sorting key function for a list of genotypes according to desired
        genotype retention priority
        """

        g = gt.get('GT', (0, ))
        is_nocall = None in g or '.' in g
        g = tuple([int(a) for a in g if a is not None and a != '.'])
        ac = sum(g)
        is_nonref = ac > 0
        tb = gt.get(tiebreak, -10e10)
        if tb is None:
            tb = -10e10

        return is_nocall, is_nonref, tb

    # Integrate each sample's genotypes across records
    for sid in all_sids:
        gts = [r.samples[sid] for r in records]
        newgts[sid].clear()
        newgts[sid].update(sorted(gts, key=__sort_gts, reverse=True)[0])

    return target_record


def integrate_infos(records):
    """
    Integrate the INFO fields of two or more pysam.VariantRecords
    Takes the union for all non-numeric fields
    Tries to be intelligent about handling numeric values based on key name (e.g., MIN, MAX)
    Takes mean of numeric values if best behavior is not obvious
    Skips protected values like END, AC, AN, AF
    Implements specialized logic for GATK-SV's CPX_INTERVALS
    Returns : pysam.libcbcf.VariantRecordInfo
    """

    if len(records) == 1:
        return records[0].info

    # Arbitrarily initialize new info as a cleared copy of the first record
    rtemp = records[0].copy()
    newinfo = rtemp.info
    newinfo.clear()

    # Get all keys present in any records
    keys = set(recursive_flatten([r.info.keys() for r in records]))

    # Exclude protected keys
    protected_keys = 'END AC AN AF'.split()
    for pk in protected_keys:
        keys = [k for k in keys if not k.startswith(pk) and not k.endswith(pk)]

    def __resolve_numerics(key, vals, key_type='Float'):
        """
        Helper function to handle smart resolution of numeric values
        """
        if key.upper().startswith('MIN') or key.upper().endswith('MIN'):
            res = np.nanmin(vals)
        if key.upper().startswith('MAX') or key.upper().endswith('MAX'):
            res = np.nanmax(vals)
        else:
            res = np.nanmean(vals)
        if key_type == 'Integer':
            return int(res)
        else:
            return float(res)

    # Iterate over all keys and update
    for key in keys:

        # Get all values to integrate
        vals = [v for v in [r.info.get(key) for r in records] if v is not None]

        # Determine handling of each key based on its corresponding header record
        key_h = [h for h in [r.header.info.get(key) for r in records] if h is not None][0]
        key_is_numeric = key_h.type in 'Integer Float'.split()
        key_is_flag = key_h.type == 'Flag'
        key_is_iterable = any([isinstance(v, Iterable) and not isinstance(v, str) for v in vals])

        # Handle integration based on value length and type
        if key_is_iterable:

            # Check for fixed length iterables if necessary
            all_same_len = len(set(map(len, vals))) == 1
            if not all_same_len:
                if key_h.number != '.':
                    msg = 'Key {} is reported as Number={} in the header, but ' + \
                          'the INFO values vary in length across records to be ' + \
                          'integrated ({})'
                    stop(msg.format(key, key_h.number, ', '.join([r.id for r in records])))
                elif key_is_numeric:
                    msg = 'Key {} is reported as Type={} in the header, but ' + \
                          'variable lengths of INFO values make integration of ' + \
                          'this field to ambiguous to proceed for these records ({})'
                    stop(msg.format(key, key_h.type, ', '.join([r.id for r in records])))

            out_type = multimode([type(v) for v in vals])[0]
            
            # Special handling of CPX_INTERVALS
            if key == 'CPX_INTERVALS':
                newinfo[key] = \
                    integrate_cpx_intervals(records,
                                            records[0].get('CPX_TYPE'),
                                            int(np.floor(np.nanmedian([r.pos for r in records]))),
                                            int(np.floor(np.nanmedian([r.stop for r in records]))))

            # Otherwise, handle integration of all other list-like/iterable INFO values
            elif key_is_numeric:
                nvl = []
                for i in range(len(vals[0])):
                    nvl.append(__resolve_numerics(key, [v[0] for v in vals], key_h.type))
                newinfo[key] = out_type(nvl)
            else:
                newinfo[key] = out_type(sorted(list(set(recursive_flatten(vals)))))

        else:
            if key_is_numeric:
                newinfo[key] = __resolve_numerics(key, vals, key_h.type)

            elif key_is_flag:
                if any(vals):
                    newinfo[key] = True

            else:
                newinfo[key] = sorted(multimode(vals))[0]

    return newinfo


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


def name_record(record, suffix_length=10):
    """
    Wrapper for name_variant designed for pysam.VariantRecord input
    """

    ref, alt = record.alleles[:2]
    if 'SVLEN' in record.info.keys():
        varlen = int(record.info['SVLEN'])
    else:
        varlen = np.abs(len(alt) - len(ref))
    vc, vsc = classify_variant(ref, alt, varlen)
    return name_variant(record.chrom, record.pos, ref, alt, vc, vsc, varlen)


def is_multiallelic(record):
    """
    Check if pysam.VariantRecord is multiallelic (including mCNVs)
    """

    if 'MULTIALLELIC' in record.filter \
    or len(record.alleles) > 2 \
    or record.info.get('SVTYPE', '') in 'CNV MCNV'.split():
        return True
    else:
        return False


