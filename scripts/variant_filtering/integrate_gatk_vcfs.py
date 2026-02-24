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
import g2cpy
import gzip
import numpy as np
import pysam


vcs = 'indel sv'.split()
newflag = '##INFO=<ID=INTEGRATED_INDEL_SV,Number=0,Type=Flag,Description=' + \
          '"Variant represents the integration of a large indel from GATK-HC ' + \
          'and small SV from GATK-SV">'


def next_record(vcf):
    """
    Helper wrapper to draw the next record from a VCF or else return None
    """
    try:
        return next(vcf)
    except:
        return None


def inputs_exhausted(pointers):
    """
    Returns true if both input pointers are exhausted
    """

    return all([r is None for r in pointers.values()])


def next_vc(pointers):
    """
    Returns next vc to process based on coordinate
    """

    pos = {vc : record.pos for vc, record in pointers.items() if record is not None}
    return sorted(pos.items(), key=lambda x: x[1])[0][0]


def integrate_records(indels, svs, header):
    """
    Integrate indel and SV records
    """

    n_indels = len(indels)
    n_svs = len(svs)

    # We always assume GATK-HC indel coordinates are more precise
    # We take the median start & end as the final record
    start = int(np.floor(np.nanmedian([r.pos for r in indels])))
    end = int(np.ceil(np.nanmedian([r.stop for r in indels])))

    # Determine which record should be used as the base
    # If multiple indels match the same SV, merge them all into that SV for simplicity
    # Otherwise, we assume the indel is more accurate
    if n_indels > 1 and n_svs == 1:
        new_rec = svs[0].copy()
    else:
        new_rec = indels[0].copy()
    new_rec.translate(header)
    new_rec.start = start
    new_rec.stop = end

    # Retain the union of all FILTERs
    for other in indels + svs:
        for f in other.filter.keys():
            new_rec.filter.add(f)

    # Update INFO
    vc, vsc, varlen = g2cpy.classify_record(new_rec, return_varlen=True)
    new_info = g2cpy.integrate_infos(indels + svs, header=header)
    new_rec.info.clear()
    new_rec.info.update(new_info)
    new_rec.info['INTEGRATED_INDEL_SV'] = True

    # Selectively update certain features depending on final variant class
    new_rec.info['ALGORITHMS'] += ('gatkhc', )
    if vc == 'sv':
        new_rec.info['SVLEN'] = int(varlen)
        new_rec.info['SVTYPE'] = vsc.upper()
        # We treat all GATK-HC variants as SR evidence due to GATK-HC local realignments
        if 'SR' not in new_rec.info.get('EVIDENCE', tuple()):
            new_rec.info['EVIDENCE'] += ('SR', )
    else:
        for k in 'SVLEN SVTYPE END ALGORITHMS EVIDENCE':
            new_rec.info.pop(k)

    # Integrate genotypes
    new_rec = g2cpy.integrate_gts(new_rec, indels + svs, 
                                  nocalls_first=False, ref_second=False,
                                  header=header)

    # Update frequency annotations
    af_stats = g2cpy.compute_allele_freq_stats(new_rec)
    new_rec.info['AC'] = int(af_stats['AC'])
    new_rec.info['AN'] = int(af_stats['AN'])
    new_rec.info['AF'] = float(af_stats['AF'])

    # Assign new record ID
    new_rec.id = g2cpy.name_record(new_rec)

    return new_rec, vc


def add_sv_fields_to_indel(record, vsc, varlen):
    """
    Adds GATK-SV expected INFO annotations to GATK-HC indel records
    """

    record.stop = int(record.pos + varlen)
    record.info['SVTYPE'] = vsc.upper()
    record.info['SVLEN'] = int(varlen)
    record.info['ALGORITHMS'] = ('gatkhc', )
    record.info['EVIDENCE'] = ('SR', )

    return record


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--indel-vcf', metavar='VCF', required=True,
                        help='GATK-HC indel VCF. Required.')
    parser.add_argument('--sv-vcf', metavar='VCF', required=True,
                        help='GATK-SV indel VCF. Required.')
    parser.add_argument('--clusters', metavar='.tsv', 
                        help='Two-column .tsv of SVs and indels to be merged.')
    parser.add_argument('--out-vcf-header', metavar='VCF', required=True,
                        help='Header-only VCF for output VCFs. Required.')
    parser.add_argument('--out-prefix', metavar='path', default='./integrated',
                        help='Path and file prefix for output VCFs.')
    args = parser.parse_args()

    # Load header for output VCFs
    header = pysam.VariantFile(args.out_vcf_header).header
    header.add_line(newflag)

    # Open connections to indel and SV VCFs
    in_vcfs = {'indel' : pysam.VariantFile(args.indel_vcf),
               'sv' : pysam.VariantFile(args.sv_vcf)}

    # Read clusters to be merged
    clusters = {}
    vid_lookups = {vc : dict() for vc in vcs}
    if args.clusters is not None:
        if 'compressed' in g2cpy.determine_filetype(args.clusters):
            cfin = gzip.open(args.clusters, 'rt')
        else:
            cfin = opne(args.clusters, 'rt')
        for i, line in enumerate(cfin.readlines()):
            if line.startswith('#'):
                continue
            sv, indel = line.rstrip().split('\t')
            sv_ids = set(sv.split(','))
            indel_ids = set(indel.split(','))
            clusters[i] = {'indel' : {'vids': indel_ids, 
                                      'n' : len(indel_ids),
                                      'records' : []},
                           'sv' : {'vids': indel_ids, 
                                   'n' : len(indel_ids),
                                   'records' : []},}
            for vid in sv_ids:
                vid_lookups['sv'][vid] = i
            for vid in indel_ids:
                vid_lookups['indel'][vid] = i
        cfin.close()

    # Load first record from indel and SV VCFs
    pointers = {vc : next_record(in_vcfs[vc]) for vc in vcs}

    # Open connections to output VCFs
    out_vcfs = {vc : pysam.VariantFile('{}.{}.vcf.gz'.format(args.out_prefix, vc), 
                                       'w', header=header) for vc in vcs}

    # Iterate until both input VCFs are exhausted
    while not inputs_exhausted(pointers):
        
        # Determine next record to process
        cur_vc = next_vc(pointers)
        cur_rec = pointers[cur_vc]

        # Check if record should be held in memory to be collapsed
        if cur_rec.id in vid_lookups[cur_vc].keys():
            i = vid_lookups[cur_vc][cur_rec.id]
            clusters[i][cur_vc]['records'].append(cur_rec)

            # Check if this record was the final one necessary to merge and output cluster
            if clusters[i]['indel']['n'] == len(clusters[i]['indel']['records']) \
            and clusters[i]['sv']['n'] == len(clusters[i]['sv']['records']):
                
                # Integrate records and route to appropriate output VCF
                new_rec, out_vc = integrate_records(clusters[i]['indel']['records'],
                                                    clusters[i]['sv']['records'],
                                                    header)
                if new_rec.info['AC'][0] > 0:
                    out_vcfs[out_vc].write(new_rec)

                # Clear this cluster to relieve memory pressure
                clusters.pop(i)

        # Otherwise, determine output variant class and route accordingly
        else:
            out_vc, vsc, varlen = g2cpy.classify_record(cur_rec, return_varlen=True)

            # Need to create new record against output header
            new_rec = cur_rec.copy()
            new_rec.translate(header)

            # Reformat indels being reclassified as SVs
            if out_vc == 'sv' and cur_vc == 'indel':
                new_rec = add_sv_fields_to_indel(new_rec, vsc, varlen)

            # Route new record to appropriate output VCF
            out_vcfs[out_vc].write(new_rec)

        # Update pointer
        pointers[cur_vc] = next_record(in_vcfs[cur_vc])

    # Close connections to output VCFs
    for outfile in out_vcfs.values():
        outfile.close()

if __name__ == '__main__':
    main()

