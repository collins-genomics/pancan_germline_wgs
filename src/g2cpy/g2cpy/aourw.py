#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Helper functions for data processing on the All of Us Researcher Workbench
"""


import subprocess
from re import sub
from sys import stdout, stderr


def abort_workflow(workflow_id, max_retries=20, timeout=30):
    """
    Abort a Cromwell workflow
    Returns a boolean indicating whether the abort was successful
    """

    crom_query_res = ''
    attempts = 0
    cmd = 'cromshell --no_turtle -t ' + str(timeout) + ' -mc abort ' + workflow_id
    while attempts < max_retries:
        try:
            res = subprocess.run(cmd, shell=True, check=True, text=True,
                                 capture_output=True)
            return True
        except:
            attempts += 1
    
    if attempts == max_retries:
        msg = 'Failed to abort workflow {} after {} retries\n'
        stderr.write(msg.format(workflow_id, attempts))
        return False


def check_workflow_status(workflow_id, max_retries=20, timeout=30):
    """
    Ping Cromwell server to check status of a single workflow
    """

    crom_query_res = ''
    attempts = 0
    while attempts < max_retries:
        crom_query = 'cromshell --no_turtle -t ' + str(timeout) + ' -mc status ' + workflow_id
        crom_query_res = subprocess.run(crom_query, capture_output=True, 
                                        shell=True, check=False, text=True).stdout
        try:
            res_splits = [sub('[,"]', '', s.split('":"')[1]) 
                          for s in crom_query_res.split('\n') 
                          if 'status' in s]
            if len(res_splits) > 0:
                return res_splits[0].lower()
            else:
                attempts += 1
        except:
            attempts += 1
    
    if attempts == max_retries:
        msg = 'Failed to get workflow status for {} after {} retries\n'
        stderr.write(msg.format(workflow_id, attempts))
        return 'status_check_failed'


def check_gcp_uris(uris, return_uri_list=False):
    """
    Check if a list of Google URIs exists
    If return_uri_list is True, will return a list of all URIs found
    Otherwise, will return a simple True/False if *all* URIs were found
    """

    # Build URI query
    query = 'gsutil -m ls ' + ' '.join(uris)
    query_res = subprocess.run(query, capture_output=True, shell=True, 
                               check=False, text=True)
    uris_found = query_res.stdout.rstrip().split('\n')

    # Return False unless every URI is found
    if return_uri_list:
        return uris_found
    else:
        return len(set(uris).intersection(set(uris_found))) == len(uris)


def collect_gcp_garbage(uris, dumpster_path='uris_to_delete.list'):
    """
    Writes a list of URIs to delete to a dumpster (flat text) file
    Note that the input URIs need not exist, and can include wildcards/double wildcards
    """

    garbage_raw = subprocess.run('gsutil -m ls ' + ' '.join(uris), 
                                 check=False, capture_output=True, 
                                 shell=True, text=True).stdout
    with open(dumpster_path, 'a') as fout:
        for uri in garbage_raw.rstrip().split('\n'):
            if uri.startswith('gs://'):
                fout.write(uri + '\n')


def collect_workflow_trash(workflow_ids, bucket, wdl_name, dumpster_path):
    """
    Add all execution files associated with execution of workflow_ids to the list
    of files to be deleted (the "dumpster")
    """

    # Write gsutil-compliant search strings for identifying files to delete
    oldest_ex_fmt = '{}/cromwell/execution/{}/{}/**'
    oldest_ex_uris = [old_ex_fmt.format(bucket, wdl_name, wid) for wid in workflow_ids]
    old_ex_fmt = '{}/cromwell-execution/{}/{}/**'
    old_ex_uris = [new_ex_fmt.format(bucket, wdl_name, wid) for wid in workflow_ids]
    new_ex_fmt = '{}/workflows/cromwell-executions/{}/{}/**'
    new_ex_uris = [new_ex_fmt.format(bucket, wdl_name, wid) for wid in workflow_ids]
    ex_uris = oldest_ex_uris + old_ex_uris
    old_out_fmt = '{}/cromwell/outputs/{}/{}/**'
    old_out_uris = [old_out_fmt.format(bucket, wdl_name, wid) for wid in workflow_ids]
    new_out_fmt = '{}/cromwell-outputs/{}/{}/**'
    new_out_uris = [new_out_fmt.format(bucket, wdl_name, wid) for wid in workflow_ids]
    out_uris = old_out_uris + new_out_uris
    all_uris = ex_uris + out_uris

    # Find list of all files present in execution or output buckets
    # and mark those files for deletion by writing their URIs to the dumpster
    collect_gcp_garbage(all_uris, dumpster_path)


def delete_uris(uris, rm_args='', verbose=False):
    """
    Deletes a list of GCP objects (or any string interpretable by gsutil rm)
    """

    if verbose:
        msg = 'Removing {}\n'
        stdout.write(msg.format(', '.join(uris)))
    subprocess.run(' '.join(['gsutil -m rm', rm_args] + uris), 
                   shell=True, text=True)


def count_vms(n_background=2):
    """
    Count number of active GCP VMs (useful as a proxy for Cromwell server load)
    This will also subtract n_background from the total count to reflect baseline
    VMs used by Cromwell or AoU RW all the time and unrelated to user-generated
    server load. Based on empirical experience, this number appeared to be 4
    when using AOU RW v1.0 default configurations and now appears to be 2 in
    the new Verily Pre RW v2.0 environment
    """

    gcheck = subprocess.run('gcloud compute instances list | wc -l', 
                            shell=True, capture_output=True, text=True,
                            encoding='utf-8')
    return int(gcheck.stdout) - n_background


def relocate_uri(src_uri, dest_uri, action='cp', check_exit_codes=True, verbose=False):
    """
    Moves or copies a GCP object from a source URI to a destination URI
    """

    if action not in 'cp mv'.split():
        msg = 'Action "{}" is not recognized. Please specify either "cp" or "mv".'
        exit(msg.format(action))

    if verbose:
        msg = 'Relocating {} to {}\n'
        stdout.write(msg.format(src_uri, dest_uri))
    subprocess.run(' '.join(['gsutil -m', action, src_uri, dest_uri]), 
                   shell=True, check=check_exit_codes)

