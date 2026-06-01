#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Filter and disambiguate candidate indel/SV overlaps
This is a subroutine called within the larger indel/SV integration workflow
"""


import argparse
import math
import networkx as nx
import pandas as pd
import re
from sys import stdout


def _calc_weight(pr):
    """
    Compute edge weight
    """

    return math.dist(tuple(pr['dist jaccard'.split()].values), (0, 1, ))


def count_prefix(nodes, prefix):
    """
    Count number of nodes with a certain string prefix
    """

    return sum(1 for n in nodes if n.startswith(prefix))


def count_vcs(cluster):
    """
    Count variant class memberships for a connected cluster
    Returns: tuple of (# indels, # SVs)
    """

    nodes = list(cluster)
    n_sv = count_prefix(nodes, "sv_")
    n_indel = count_prefix(nodes, "indel_")

    return tuple([n_indel, n_sv])


def cluster_is_valid(cluster):
    """
    Checks whether a cluster of SVs/indels is valid (is singular for SVs or indels)
    """

    n_indel, n_sv = count_vcs(cluster)

    return (n_indel == 1) or (n_sv == 1)


def prune_cluster(G, cluster):
    """
    Remove the weakest edge (largest cdist) from a given cluster
    """

    subG = G.subgraph(cluster).copy()
    
    # Sort edges by descending weight (larger weight = worse match)
    edges_sorted = sorted(
        subG.edges(data=True),
        key=lambda x: x[2].get("weight", 0),
        reverse=True
    )

    # Remove the single worst edge
    u, v, data = edges_sorted[0]
    G.remove_edge(u, v)


def iterative_prune(G):
    """
    Prunes all clusters in an entire indel/SV graph
    """

    while True:

        removed_total = 0
        
        # Recompute connected components each round
        clusters = list(nx.connected_components(G))
        
        for cluster in clusters:
            if not cluster_is_valid(cluster):
                prune_cluster(G, cluster)
                removed_total += 1
        
        if removed_total == 0:
            break
    
    return G


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--pair-distance', metavar='.tsv', required=True,
                        help='Three-column .tsv of variant IDs and distance per pair')
    parser.add_argument('--pair-jaccard', metavar='.tsv', required=True,
                        help='Three-column .tsv of variant IDs and genotype jaccard per pair')
    parser.add_argument('--minimum-jaccard', type=float, default=0,
                        help='Minimum Jaccard index [default: 0]', 
                        metavar='[float]')
    parser.add_argument('-o', '--outfile', metavar='[path]', default='stdout',
                        help='Path to output .tsv of final variant clusters')
    args = parser.parse_args()

    # Load pairs and Euclidean distances
    pairs = pd.read_csv(args.pair_distance, sep='\t', header=None,
                        names='vid1 vid2 dist'.split())

    # Load Jaccards, filter, and join to distances
    jac = pd.read_csv(args.pair_jaccard, sep='\t', header=0)
    jac.columns = jac.columns.str.replace('#', '')
    jac = jac[jac.jaccard >= args.minimum_jaccard]
    pairs = pd.merge(pairs, jac, how='inner')

    # Compute Euclidean distance between all pairs when accounting for Jaccard index
    pairs['weight'] = pairs.apply(_calc_weight, axis=1)

    # Build graph of all pairs
    G = nx.from_pandas_edgelist(pairs, source="vid1", target="vid2", edge_attr="weight")

    # Iterate over all subgraphs until no subgraph contains multiple indels *and* SVs    
    G = iterative_prune(G)

    # Open connection to output file and write header
    if args.outfile in '- stdout /dev/stdout'.split():
        fout = stdout
    else:
        fout = open(args.outfile, 'w')
    fout.write('#SVs\tindels\n')

    # Write all pruned clusters to the output file
    for cluster in nx.connected_components(G):
        vids = list(cluster)
        svs = [re.sub('^sv_', '', vid) for vid in vids if vid.startswith("sv_")]
        indels = [re.sub('^indel_', '', vid) for vid in vids if vid.startswith("indel_")]
        if len(svs) > 0 and len(indels) > 0:
            fout.write('{}\t{}\n'.format(','.join(sorted(svs)), ','.join(sorted(indels))))

    # Close outfile to clear buffer
    fout.close()


if __name__ == '__main__':
    main()

