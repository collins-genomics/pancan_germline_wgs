#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
VCF parsing & manipulation functions
"""


import math
import numpy as np
import re
from collections import Counter
from collections.abc import Iterable
from functools import partial
from statistics import multimode
from .genomics import classify_variant, name_variant
from .utilities import recursive_flatten


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


def classify_record(record, return_varlen=False):
    """
    Wrapper for classify_variant designed for pysam.VariantRecord input
    """

    ref, alt = record.alleles[:2]
    if 'SVLEN' in record.info.keys():
        varlen = int(record.info['SVLEN'])
    else:
        varlen = np.abs(len(alt) - len(ref))
    outs = list(classify_variant(ref, alt, varlen))
    if return_varlen:
        outs.append(varlen)
    return tuple(outs)


def compute_allele_freq_stats(record, keys=['AC', 'AN', 'AF', 'N_REF', 'N_HET', 'N_HOM', 'N_MISSING']):
    """
    Helper wrapper around apply_across_samples() to compute common frequency 
    stats for a pysam.VariantRecord object
    Returns: dict of floats keyed by 'keys'
    """

    res = {k : 0 for k in keys}

    for sid in record.samples:
        # Parse GT field
        av = parse_gt(record.samples[sid]['GT'])
        ac = int(av.get('AC', 0))
        an = int(av.get('AN', 0))
        mis = int(av.get('N_missing', 0))
        ref = an - mis - ac

        # Increment allele counts
        res['AC'] += ac
        res['AN'] += int(np.max([an - mis, 0]))

        # Parse genotype
        if ac == 0 and ref > 0:
            res['N_REF'] += 1
        elif ac == an:
            res['N_HOM'] += 1
        elif ac > 0 and ref > 0:
            res['N_HET'] += 1
        elif mis == an:
            res['N_MISSING'] += 1
        else:
            print('Unable to parse GT {}\n'.format(record.samples[sid]['GT']))
            import pdb; pdb.set_trace()

    res['AF'] = float(res['AC'] / res['AN'])

    all_keys = list(res.keys())
    for ak in all_keys:
        if ak not in keys:
            res.pop(ak)

    return res


def convert_gt(gt):
    """
    Converts a VCF GT string into pysam-style tuple representation, or vice versa
    """

    if isinstance(gt, tuple):
        res = []
        for a in gt:
            if a is None:
                res.append('.')
            else:
                res.append(str(a))
        conv = '/'.join(res)

    elif isinstance(gt, str):
        res = []
        for a in re.split('[/|\\|]', gt):
            if a == '.':
                res.append(None)
            else:
                res.append(int(a))
        conv = tuple(sorted(res, key=lambda x: (x, x is None)))

    return conv


def integrate_gts(target_record, records, consensus=True, sort_key='GQ', 
                  nocalls_first=True, ref_second=True,
                  nonref_override_metric=10e10, ref_override_metric=10e10,
                  collapse_formats=False, header=None):
    """
    Merge genotypes across pysam.VariantRecords

    By default, GT selection per sample is prioritized as:
    1. Modal GT among `records` (only if records is a nested list and consensus=True, otherwise ignored)
    2. Null / no-call (./.)
    3. All non-ref GTs (moved ahead of null/no-call if non-ref sort_key 
       is greater than nonref_override_metric)
    4. All ref GTs (moved ahead of null/no-call if ref sort_key is greater than ref_override_metric)
       Ties are broken by sort_key, and entire FORMAT is retained from winner GT 
       (unless collapse_formats is True, in which case all records with the winning 
        GT will have their info added overwritten in ascending order by sort_key)

    Writes new genotypes into target_record and returns updated record
    """

    nested_input = any(isinstance(item, list) for item in records)

    if len(records) < 2 and not nested_input:
        return target_record

    # Reheader all records, if necessary
    if not nested_input:
        records = [records]
    if header is not None:
        for i in range(len(records)):
            records[i] = [r.copy() for r in records[i]]
            for r in records[i]:
                r.translate(header)
        fmt_meta = {f: header.formats[f] for f in header.formats.keys()}
    else:
        fmt_meta = {f: target_record.header.formats[f] for f in target_record.format.keys()}

    # Check to ensure all samples are the same across all records
    flat_records = recursive_flatten(records)
    sids_per_rec = [set(r.samples.keys()) for r in flat_records]
    all_sids = set(recursive_flatten(sids_per_rec))
    if len(set(map(len, [all_sids] + sids_per_rec))) > 1:
        msg = 'Attempted to integrate genotypes but failed due to inconsistent ' + \
              'samples for records {}'
        exit(msg.format(', '.join([r.id for r in flat_records])))

    # Direct access to target GT fields
    newgts = target_record.samples

    def __sort_gts(gt, tiebreak='GQ', nocalls_first=nocalls_first, ref_second=ref_second,
                   ref_override=ref_override_metric, nonref_override=nonref_override_metric,
                   any_nocalls_present=False):
        """
        Custom sorting key function for a list of genotypes according to desired
        genotype retention priority
        """

        # Gather all information required to rank priority tier
        g = gt.get('GT', (0, ))
        is_nocall = None in g or '.' in g
        g = tuple([int(a) for a in g if a is not None and a != '.'])
        try:
            ac = sum(g)
        except:
            ac = 0
        tb = gt.get(tiebreak)
        if isinstance(tb, tuple):
            tb = tb[0]
        if tb is None:
            tb = -10e10
        is_ref = len(g) > 0 and ac == 0
        is_nonref = ac > 0

        # Base tier assignment
        if is_nocall:
            if nocalls_first:
                tier = 2
            else:
                tier = -1
        else:
            if is_nonref:
                tier = 1
            else:
                if ref_second:
                    tier = 0
                else:
                    tier = 1

        # Increment tiers for defined GTs if they exceed override values, but 
        # only if any no-call GTs are present
        if any_nocalls_present:
            if is_nonref and tb >= nonref_override:
                tier += 1
            if is_ref and tb >= ref_override:
                tier += 1

        return (tier, tb)

    # Integrate each sample's genotypes across records
    for sid in all_sids:

        # First, extract _all_ FORMAT entries for all records
        if nested_input:
            sdat = []
            for subrecs in records:
                sdat.append([r.samples[sid] for r in subrecs])
        else:
            sdat = [r.samples[sid] for r in records]

        # Determine which GT to keep. Behavior depends on nested_input and consensus
        if nested_input:
            if consensus:

                # Transpose nested FORMATs for easy consensus identification
                keyed_sdat = []
                for i in range(len(sdat)):
                    keyed_sdat.append(dict())
                    for sd in sdat[i]:
                        gt = sd.get('GT', (None, None,))
                        gv = parse_gt(gt).values()
                        gk = '_'.join([str(v) for v in gv])
                        if gk not in keyed_sdat[i].keys():
                            keyed_sdat[i][gk] = []
                        keyed_sdat[i][gk].append(sd)

                # Consensus GT is only defined if there is a unique modal GT (phase-agnostic)
                gt_k = Counter(recursive_flatten([list(sl.keys()) for sl in keyed_sdat]))
                best_k = max(gt_k.values())
                best_k_gts = [gt for gt, k in gt_k.items() if k == best_k]
                consensus_reached = len(best_k_gts) == 1
                # If consensus is found, subset all sdat to only those GTs
                if consensus_reached:
                    consensus_gt = best_k_gts[0]
                    sdat = [sd.get(consensus_gt, []) for sd in keyed_sdat]

            # Always flatten sdat before proceeding
            # Note: do not use recursive_flatten() here as it will extract the
            # values from the VariantRecordSample objects (since they are iterable)
            sdat = [sd for sl in sdat for sd in sl]

        any_nocalls = any([None in gt.get('GT', (None,)) for gt in sdat])
        _sorter = partial(__sort_gts, any_nocalls_present=any_nocalls)
        keep_gts = sorted(sdat, key=_sorter, reverse=True)
        if consensus and collapse_formats:
            if consensus_reached:
                keep_gts.reverse()
            else:
                keep_gts = [keep_gts[0]]
        else:
            keep_gts = [keep_gts[0]]
        
        # Next, clear existing GT and replace all values with None
        # For reasons that aren't obvious, in very rare edge cases pysam fails to
        # write None values; given this, it's weirdly safer to first try to 
        # assign values as the string '.' to force pysam to write out the VCF 
        # spec indication for a missing FORMAT value
        newgts[sid].clear()
        for dfk in target_record.format.keys():
            dkn = fmt_meta[dfk].number
            for dv in ['.', None]:
                if isinstance(dkn, int):
                    try:
                        if dkn == 1:
                            def_val = dv
                        else:
                            def_val = tuple([dv] * dkn)
                        newgts[sid][dfk] = def_val
                        break
                    except:
                        pass
                elif dkn == '.':
                    try:
                        newgts[sid][dfk] = dv
                        break
                    except:
                        pass
                else:
                    for i in range(1, 10):
                        try:
                            newgts[sid][dfk] = tuple([dv] * i)
                            break
                        except:
                            pass
        
        # Finally, update cleared GT with GT to keep
        for sd in keep_gts:
            for fk, fv in sd.items():
                try:
                    newgts[sid][fk] = fv
                except:
                    pass

    return target_record


def integrate_infos(records, do_cpx_intervals=False, header=None):
    """
    Integrate the INFO fields of two or more pysam.VariantRecords
    Takes the union for all non-numeric fields
    Tries to be intelligent about handling numeric values based on key name (e.g., MIN, MAX)
    Takes mean of numeric values if best behavior is not obvious
    Skips protected values like END, AC, AN, AF
    Implements specialized logic for GATK-SV's CPX_INTERVALS, but only if do_cpx_intervals is True
    Returns : pysam.libcbcf.VariantRecordInfo
    """

    if len(records) == 1:
        return records[0].info

    if header is not None:
        records = [r.copy() for r in records]
        for r in records:
            r.translate(header)

    # Arbitrarily initialize new info as a cleared copy of the first record
    rtemp = records[0].copy()
    newinfo = rtemp.info
    newinfo.clear()

    # Get all keys present in any records
    keys = set(recursive_flatten([r.info.keys() for r in records]))

    # Exclude protected keys
    protected_keys = 'END AC AN AF SVLEN'.split()
    for pk in protected_keys:
        keys = [k for k in keys if not k.startswith(pk) and not k.endswith(pk)]

    def __resolve_numerics(key, vals, key_type='Float'):
        """
        Helper function to handle smart resolution of numeric values
        """
        vals = [v for v in vals if v is not None]
        if len(vals) == 0:
            return None
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
            if key == 'CPX_INTERVALS' and do_cpx_intervals:
                newinfo[key] = \
                    integrate_cpx_intervals(records,
                                            records[0].info.get('CPX_TYPE'),
                                            records[0].chrom,
                                            int(np.floor(np.nanmedian([r.pos for r in records]))),
                                            int(np.floor(np.nanmedian([r.stop for r in records]))))

            # Otherwise, handle integration of all other list-like/iterable INFO values
            elif key_is_numeric:
                nvl = []
                for i in range(len(vals[0])):
                    nvl.append(__resolve_numerics(key, [v[i] for v in vals], key_h.type))
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


def is_multiallelic(record):
    """
    Check if pysam.VariantRecord is multiallelic (including mCNVs)
    """

    if 'MULTIALLELIC' in record.filter \
    or len(record.alleles) > 2:
        return True
    if 'SVTYPE' in record.info.keys():
        return record.info.get('SVTYPE', '') in 'CNV MCNV'.split()
    return False


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


def parse_gt(gt):
    """
    Digests a pysam GT-style tuple into a dict keyed by AC, AN, and N_missing
    Note that AN here reflects total ploidy, not number of defined/called alleles
    This is slightly different behavior than standard VCF spec but can be
      reconciled by comparing N_missing to AN
    """

    an = len(gt)
    gt_def = [a for a in gt if a is not None]
    n_missing = an - len(gt_def)
    ac = len([a for a in gt_def if a > 0])

    return {'AC' : ac, 'AN' : an, 'N_missing' : n_missing}


