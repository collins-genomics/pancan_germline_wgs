#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Gather variant-level annotations for GT filtering model
Input VCF must strictly contain only a single variant type (snv, indel, or SVs)
"""


import argparse
import gzip
import g2cpy
import pybedtools as pbt
import pyBigWig as pbw
import pysam
from sys import stdin, stdout


all_vcs = 'snv indel sv'.split()
all_vscs = {'snv' : 'ti tv'.split(),
            'indel' : 'ins del'.split(),
            'sv' : 'DEL DUP CNV INS CPX OTH'.split()}


def gather_site_features(record, n_samples, filters, fbts=dict(), 
                         fbws=dict(), colnames=False):
    """
    Main function to gather and format site features and, optional, feature names
    """

    outcols = None

    # Get basic descriptives
    vid = record.id
    qual = float(record.qual)
    outvals = [vid, qual]
    if colnames:
        outcols = 'vid qual'.split()

    # One-hot encode variant subclass
    # For simplicity, recode all inversions as CPX for purposes of filtering
    # And recode all unusual SVs into a catch-all OTH category
    vc, vsc, varlen = g2cpy.classify_record(record, return_varlen=True)
    if vsc == 'INV':
        vsc = 'CPX'
    if vc == 'sv' \
    and ( vsc in 'CNV mCNV MCNV'.split() \
          or 'MULTIALLELIC' in record.filter.values() ):
        vsc = 'CNV'
    if vc == 'sv' and vsc not in all_vscs[vc]:
        vsc = 'OTH'
    for k in all_vscs[vc]:
        if vsc == k:
            outvals.append(1)
        else:
            outvals.append(0)
    if colnames:
        outcols += all_vscs[vc]
    if vc == 'snv':
        varlen = None
    else:
        outvals.append(int(varlen))
        if colnames:
            outcols.append('varlen')

    # Annotate non-PASS FILTERs
    for f in filters:
        if f in record.filter.values():
            outvals.append(1)
        else:
            outvals.append(0)
    if colnames:
        outcols += [k.lower() for k in filters]

    # Get variant frequency info
    af_stats = g2cpy.compute_allele_freq_stats(record)
    ac = int(af_stats.get('AC', 0))
    an = int(af_stats.get('AN', 0))
    af = float(af_stats.get('AF', 0))
    ncr = float(af_stats.get('N_MISSING', 0) / n_samples)
    hwe_stats = list(g2cpy.hwe_dist(af_stats))
    outvals += [ac, af, ncr] + hwe_stats
    if colnames:        
        outcols += 'ac af ncr hwe_d ex_het ex_hom'.split()

    # Get GATK-SV specific features
    if vc != 'snv':
        algs = ','.join(sorted(list(record.info.get('ALGORITHMS', tuple()))))
        ev = ','.join(sorted(list(record.info.get('EVIDENCE', tuple()))))
        outvals += [algs, ev]
        if colnames:        
            outcols += 'algs site_ev'.split()
        sv_bools = 'BOTHSIDES_SUPPORT HIGH_SR_BACKGROUND IMPUTED ' + \
                   'PESR_GT_OVERDISPERSION INTEGRATED_INDEL_SV'
        sv_bools = sv_bools.split()
        for k in sv_bools:
            if record.info.get(k, False):
                outvals.append(1)
            else:
                outvals.append(0)
        outvals.append(float(record.info.get('FRAC_IMPUTED', 0)))
        if colnames:
            outcols += [k.lower() for k in sv_bools]
            outcols.append('FRAC_IMPUTED'.lower())

    # Get GATK-HC specific features
    hc_feats = 'BaseQRankSum FS MQRankSum MQ ReadPosRankSum SOR'.split()
    outvals += [record.info.get(k, '.') for k in hc_feats]
    if colnames:        
        outcols += [k.lower() for k in hc_feats]

    # Add BED feature annotations, if provided
    if len(fbts) > 0:
        for fname, fbt in fbts.items():
            outvals.append(annotate_feature(record, fbt, slop=1))
            if colnames:
                outcols.append(fname)

    # Add bigWig feature annotations, if provided
    if len(fbws) > 0:
        for fname, fbw in fbws.items():
            bwv = list(annotate_feature(record, fbw, slop=75))
            if vc == 'snv':
                outvals.append(float(bwv[0]))
                if colnames:
                    outcols.append(fname)
            else:
                outvals += bwv
                if colnames:
                    outcols += [p + '_' + fname for p in 'min max'.split()]

    return outvals, outcols


def annotate_feature(record, feature, slop=0):
    """
    Annotate intersection between a variant `record` and a `feature` object

    If `feature` is pbt.BedTool:
        Annotate overlap between any features in `fbt` and the POS and END of `record`

        Returns a numeric value indicating whether zero, one, or both of POS/END
        overlap with `fbt` (SNVs have max of 1, indels/SVs have max of 2)

    If `feature` is pbw.bigWigFile:
        Annotate mean of `fbw` around each of the ends of `record`
        A window of ±`slop` is drawn around each of POS and END 

        Returns a sorted tuple of average values. One value is returned for SNVs,
        while two values are returned for indels/SVs
    """

    # Enumerate variant boundary coordinates
    coords = [(record.chrom, 
               max([record.pos - slop, 0]),
               record.pos + slop)]
    vc, vsc = g2cpy.classify_record(record)
    if vc != 'snv':
        coords.append((record.chrom,
                       max([record.stop - slop]),
                       record.stop + slop))

    if isinstance(feature, pbt.BedTool):
        vstr = '\n'.join(['{}\t{}\t{}\n'.format(*v) for v in coords])
        vbt = pbt.BedTool(vstr, from_string=True)
        return len(vbt.intersect(feature, u=True))

    if hasattr(feature, "isBigWig"):
        if feature.isBigWig():
            means = [feature.stats(str(v[0]), int(v[1]), int(v[2]))[0] for v in coords]
            sorted_means = sorted([n for n in means if isinstance(n, float)])
            if len(sorted_means) < len(means):
                sorted_means += ['.'] * (len(means) - len(sorted_means))
            return tuple(sorted_means)

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
                        type=str, help='Output .tsv')
    parser.add_argument('--feature-bed', metavar='feature_name=path', action='append',
                        help='BED-style track to annotate binary overlap with ' +
                        'variant boundaries (±1bp). Can be provided any number of times. ' +
                        'Must be specified as feature_name=path/to/features.bed.')
    parser.add_argument('--feature-bigwig', metavar='feature_name=path', action='append',
                        help='bigWig-style track to annotate mean at variant ' +
                        'boundaries (±75bp). Can be provided any number of times. ' +
                        'Must be specified as feature_name=path/to/feature.bw.')
    parser.add_argument('-p', '--precision', default=3, metavar='integer',
                        type=int, help='Floating point precision')
    args = parser.parse_args()

    # Open connection to input VCF
    if args.input_vcf in 'stdin /dev/stdin -'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.input_vcf)

    # Get constants to be reused during record parsing
    n_samples = len(invcf.header.samples)
    filters = [f.name for f in invcf.header.filters.values() \
               if f.name not in 'PASS MULTIALLELIC UNRESOLVED'.split() 
               and f.name != '.']

    # Open connections to input feature BEDs, if optioned
    fbts = {}
    if len(args.feature_bed) > 0:
        for fbti in args.feature_bed:
            fname = str(fbti.split('=')[0]).lower()
            fbt = pbt.BedTool('='.join(fbti.split('=')[1:]))
            fbts[fname] = fbt
    

    # Open connections to input feature bigWigs, if optioned
    fbws = {}
    if len(args.feature_bigwig) > 0:
        for fbwi in args.feature_bigwig:
            fname = str(fbwi.split('=')[0]).lower()
            fbw = pbw.open('='.join(fbwi.split('=')[1:]))
            fbws[fname] = fbw

    # Open connection to output TSV
    if args.output_tsv in 'stdin /dev/stdin -'.split():
        outfile = stdout
    else:
        if 'compressed' in g2cpy.determine_filetype(args.output_tsv):
            outfile = gzip.open(args.output_tsv, 'wt', encoding='utf-8')
        else:
            outfile = open(args.output_tsv, 'w')

    # Iterate over input VCF and collect variant annotations for each
    for i, record in enumerate(invcf):

        outvals, outcols = gather_site_features(record, n_samples, filters, 
                                                fbts, fbws, colnames=(i == 0))

        # Write header if first line
        if i == 0:
            outfile.write('\t'.join(outcols) + '\n')

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

