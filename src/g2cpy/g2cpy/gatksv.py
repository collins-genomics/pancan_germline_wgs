#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
GATK-SV utility functions
"""

import pybedtools as pbt
import re
from .utilities import recursive_flatten


def integrate_cpx_intervals(records, cpx_type, cpos, cend):
    """
    Integrates INFO/CPX_INTERVALS across two or more pysam.VariantRecord objects
    corresponding to SVTYPE=CPX from GATK-SV

    Returns tuple of integrated CPX_INTERVALS that can be assigned to a single record
    """

    itypes = 'del dup ins inv'.split()

    # First, infer how many intervals there _should_ be
    exp_ints = {itype : cpx_type.count(itype) for itype in itypes}

    # Next, count how many intervals there _actually_ are
    cpx_ints = tuple(set(recursive_flatten([r.get('CPX_INTERVALS') for r in records])))
    obs_ints = {itype : sum([x.lower().startswith(itype) for x in cpx_ints]) 
                for itype in itypes}

    # Don't worry about inverted orientations for dispersed duplications
    if 'ddup' in cpx_type or 'ins' in cpx_type:
        exp_ints.pop('inv')
        obs_ints.pop('inv')

    # If interval types match, do nothing
    if exp_ints == obs_ints:
        return cpx_ints

    # Otherwise, clean up intervals s/t there are no more than expected
    recbt = pbt.BedTool('{}\t{}\t{}\n'.format(records[0].chrom, cpos, cend), 
                        from_string=True)
    for itype in itypes:

        # Skip interval types with the same or fewer than expected intervals
        if obs_ints[itype] <= exp_ints[itype]:
            continue

        # Blanket prune unexpected interval types
        if exp_ints[itype] == 0:
            cpx_ints = [i for i in cpx_ints if not i.startswith(itype.upper())]
            continue

        # Otherwise, try merging overlapping intervals with increasingly permissive
        # overlap requirements until the expected number of intervals are retained
        ibt_str = [re.sub('[:-]', '\t', i.split('_')[1]) 
                   for i in cpx_ints if i.startswith(itype.upper())]
        ibt = pbt.BedTool('\n'.join(ibt_str), from_string=True)
        # TODO: implement this

        # OLD
        # idf = ibt.coverage(recbt).to_dataframe()
        # idf = idf.sort_values(by=idf.columns[-1], ascending=False)
        # keep_ints = ['{}_{}:{}-{}'.format(itype.upper(), *il) 
        #              for il in idf.iloc[0:exp_ints[itype], 0:3].values.tolist()]
        # cpx_ints = [i for i in cpx_ints if not i.startswith(itype.upper())] + keep_ints

    # Trim all intervals s/t they do not extend beyond cpos/cend
    # TODO: implement this

    return tuple(cpx_ints)

