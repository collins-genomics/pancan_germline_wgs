#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Routine to submit, track, manage, and stage a method parallelized per chromosome
"""


# Import libraries
import argparse
import g2cpy
import json
import numpy as np
import pandas as pd
import subprocess
import textwrap
import threading
from datetime import datetime
from os import environ, getenv, path, putenv, remove
from re import sub
from sys import stdout, stderr, exit
from time import sleep


# Define constants
hg38_primary_contigs = ['chr{}'.format(x+1) for x in range(22)] + ['chrX', 'chrY']
strict_sub_status = 'failed aborted aborting unstaged not_started'.split()
running_status = 'running doomed submitted'.split()

# Define global thread lock
print_lock = threading.Lock()


def clean_date():
    """
    Wrapper for clean date & time reporting
    """

    return datetime.now().strftime("%Y-%m-%d at %H:%M:%S")


def startup_report(inputs):
    """
    Prints a startup message to stdout with all runtime info
    """

    start_log="""

    #=======================================================#
    |      DFCI G2C chromosome-sharded workflow manager     |
    | Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu> |
    #=======================================================#

    Launched on {}

    Inputs:"""
    start_log = start_log.format(clean_date())
    print(textwrap.dedent(start_log))

    for key, val in inputs.items():
        if val is not None:
            if isinstance(val, bool):
                if val:
                    print('  - {} enabled'.format(key))
            else:
                print('  - {}: {}'.format(key, val))
    print('')


def report_status(all_status, title=None):
    """
    Summarizes workflow status across all chromosomes held in all_status tracker
    """

    # Group contigs by status
    status_groups = {}
    for k, status in all_status.items():
        if status not in status_groups:
            status_groups[status] = set()
        status_groups[status].add(k)

    # Report status in decreasing frequency
    if title is None:
        print('[{}] Status summary for all contigs:'.format(clean_date()))
    else:
        print(title)
    for status, contigs in sorted(status_groups.items(), 
                                  key=lambda x: len(x[1]), reverse=True):
        print('  - {:,} {}: {}'.format(len(contigs), status, 
                                       ', '.join(g2cpy.chromsort(list(contigs)))))
    print('')


def vm_load_daemon(period_mins=5, n_retries=10):
    """
    Background thread that reports current GCP VM load every `period_mins` minutes
    """

    msg = '[{}] Current Cromwell server load: {:,} active VMs'

    while True:
        with print_lock:
            print(msg.format(clean_date(), g2cpy.count_vms()))
        sleep(60 * period_mins)


def relocate_outputs(workflow_id, staging_bucket, wdl_name, output_json_uri,
                     action='cp', timeout=30, max_retries=20, verbose=False):
    """
    Relocate final output files for a single workflow to a permanent storage location
    """

    msg = 'Relocating {} to {}\n'

    # Use cromshell list-outputs to get and read a melted dict-like .txt of outputs
    crom_query = 'cromshell --no_turtle -t ' + str(timeout) + \
                 ' --machine_processable list-outputs --json-summary ' + \
                 workflow_id
    attempts = 0
    while attempts < max_retries:
        crom_query_res = subprocess.run(crom_query, capture_output=True, shell=True, 
                                        check=False, text=True).stdout
        if crom_query_res != '':
            break
        else:
            attempts += 1
        if attempts == max_retries:
            msg = 'cromshell list-outputs failed after {:,} retries with {:,}s timeout'
            exit(msg.format(max_retries, timeout))

    # Iterate over each output mentioned by cromwell and move each to staging_bucket
    # While relocating output files, build an outputs .json to reflect their new locations
    outputs_dict = {}
    for key, src_vals in json.loads(crom_query_res).items():

        if src_vals is None:
            continue

        if isinstance(src_vals, str):
            src_vals = [src_vals]
        elif isinstance(src_vals, list):
            src_vals = g2cpy.recursive_flatten(src_vals)
        else:
            msg = 'relocate_outputs cannot parse output "{}" due to unknown type'
            exit(msg.format(key))

        for src_uri in src_vals:
            
            # Explicit skipping of Nones is necessary when src_vals is rendered as a generator
            if src_uri is None:
                continue

            # Format destination URI to keep nested structure of subworkflows
            # but to remove all arbitrary Cromwell hashes
            dest_parts = [sub('^call-', '', x) for x in src_uri.split('/') 
                          if x.startswith('call-') or x.startswith('shard-')]
            dest_uri = '/'.join([staging_bucket] + dest_parts + [path.basename(src_uri)])
        
            # Stage output
            g2cpy.relocate_uri(src_uri, dest_uri, verbose=verbose)

            # Update outputs .json
            if key not in outputs_dict.keys():
                outputs_dict[key] = list()
            outputs_dict[key].append(dest_uri)

    # Write updated output .json to file and copy that file into staging_bucket
    for key in outputs_dict.keys():
        if len(outputs_dict[key]) == 1:
            outputs_dict[key] = outputs_dict[key][0]
    json_fout = path.basename(output_json_uri)
    with open(json_fout, 'w') as fout:
        json.dump(outputs_dict, fout)
    g2cpy.relocate_uri(json_fout, output_json_uri)
    remove(json_fout)


def submit_workflow(contig, wdl, input_template, input_json, options_json, 
                    prev_wids, dependencies_zip=None, contig_var_overrides=None, 
                    dry_run=False, quiet=False, verbose=False, validate=False):
    """
    Submit a workflow for a single contig to Cromwell and log the workflow ID
    """

    # Assign contig to all protected environment variables
    for vbl in 'CHR CHROM CONTIG chr chrom contig'.split():
        environ[vbl] = contig

    # Substitute contig-specific environment variables provided as optional .json
    if contig_var_overrides is not None:
        with open(contig_var_overrides) as fin:
            cvo = json.load(fin)
            if contig in cvo.keys():
                for key, value in cvo[contig].items():
                    # Note: need to replace all single quotes with double quotes
                    # to remain consistent with .json spec. This is especially
                    # an issue for arrays of strings or URIs.
                    environ[key] = sub('\'', '"', str(value))

    # Substitute variables into .json template
    with open(input_template) as fin:
        with open(input_json, 'w') as fout:
            for line in fin.readlines():
                fout.write(path.expandvars(line))

    # Build the workflow submission command
    cmd = 'cromshell --no_turtle -t 120 -mc submit'
    cmd += ' --options-json ' + options_json
    if dependencies_zip is not None:
        cmd += ' --dependencies-zip ' + dependencies_zip
    if not validate:
        cmd += ' --no-validation'
    cmd += ' ' + wdl + ' ' + input_json

    # Submit the workflow
    if dry_run:
        if not quiet:
            msg = '[{}] Would submit new workflow with the following ' + \
                  'command if not for --dry-run:\n{}'
            with print_lock:
                print(msg.format(clean_date(), cmd))
        status = 'dry_run_skipped'
    else:
        sub_res = subprocess.run(cmd, capture_output=True, shell=True, 
                                 check=False, text=True)
        if verbose and not quiet:
            msg = '[{}] Standard output from Cromwell submission command:\n{}'
            with print_lock:
                print(msg.format(clean_date(), sub_res.stdout))
        try:
            new_wid = json.loads(sub_res.stdout)['id']
            with open(prev_wids, 'a') as fout:
                fout.write(new_wid + '\n')
            status = 'submitted'
            if not quiet:
                msg = '[{}] Submitted {} for contig {} (workflow ID: {})'
                with print_lock:
                    print(msg.format(clean_date(), path.basename(wdl), contig, new_wid))
        except:
            status = 'submission_error'

    return status


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-w', '--wdl', help='WDL file to execute', required=True)
    parser.add_argument('-i', '--input-json-template', help='Template for input ' +
                        '.json to pass to Cromwell. All environment variables ' +
                        'present will be automatically expanded. Any of $CHR, ' +
                        '$CHROM, $CONTIG, $chr, $chrom, or $contig are reserverd ' +
                        'and will be replaced by contig name within this routine.',
                        required=True)
    parser.add_argument('-o', '--staging-bucket', help='G2C output staging ' + 
                        'bucket. Must provide terminal bucket within which ' +
                        'each chromosome will have outputs staged in separate ' +
                        'directory.')
    parser.add_argument('-D', '--dependencies-zip', help='WDL dependencies .zip')
    parser.add_argument('-O', '--cromwell-options-json', help='Cromwell options .zip',
                        default='code/refs/json/aou.cromwell_options.default.json')
    parser.add_argument('-V', '--contig-variable-overrides', help='Optional .json ' +
                        'of { $contig : { variable : value, ...} } pairs ' +
                        'that will be treated as environment variables for ' +
                        'substitution in --input-json-template, and will ' +
                        'supercede any global or local shell  environment ' +
                        'variables, if relevant.')
    parser.add_argument('-n', '--name', help='Method name or nickname. If ' +
                        'not provided, the basename of --wdl will be used.')
    parser.add_argument('--bucket', help='Root bucket [defaut: use ' +
                        '$WORKSPACE_BUCKET environment variable]')
    parser.add_argument('--contig-list', help='List of contig names to ' + 
                        'evaluate [default: all 24 primary autosomes]')
    parser.add_argument('-t', '--status-tsv', help='Two-column .tsv of contig ' +
                        'name and last known status')
    parser.add_argument('-r', '--base-directory', default='/home/jupyter',
                        help='Base directory for exectuion. Used to satisfy ' +
                        'assumptions about the locations of Cromshell log files ' + 
                        'and other metadata.')
    parser.add_argument('-b', '--behavior', default='patient', 
                        choices='patient strict indifferent'.split(),
                        help='Behavior for handling resubmission of failed ' +
                        'workflows. Both \'patient\' and \'strict\' will never ' +
                        'allow more than one active workflow per chromosome. ' +
                        '\'patient\' will allow any non-failed workflow to proceed ' +
                        'to completion before submitting a new one. \'strict\' ' +
                        'will instead abort a DOOMED workflow prior to submitting ' +
                        'a new one. \'indifferent\' will submit a new workflow ' +
                        'if the status of the active workflow is unclear, or as ' +
                        'soon as it becomes clear the active workflow will not ' +
                        'be successful, but takes no action to stop the active ' +
                        'workflow.')
    parser.add_argument('-d', '--dumpster', default='uris_to_delete.list',
                        help='Path to file collecting all GCP URIs to be deleted. ' +
                        'Note that this script will only ever append to this file, ' +
                        'not overwrite it.')
    parser.add_argument('--no-stage', action='store_true', help='Disable ' +
                        'automatic staging of outputs. Useful for WDL testing.')
    parser.add_argument('--no-cleanup', action='store_true', help='Disable ' +
                        'automatic cleanup of Cromwell execution buckets.')
    parser.add_argument('--hard-reset', action='store_true', help='\'Unstage\' ' +
                        'all files from each contig\'s staging bucket prior to ' +
                        'beginning monitor routine. This will force a rerun of ' +
                        'the workflow for all contigs and will also ignore all ' +
                        'prior workflow submissions when evaluating ' +
                        '--max-attempts')
    parser.add_argument('-m', '--max-attempts', type=int, default=3, 
                        help='Maximum number of submission attempts for any one contig')
    parser.add_argument('-g', '--outer-gate', type=float, default=20, 
                        help='Number of minutes to wait between monitor cycles')
    parser.add_argument('--vm-gate', type=int, default=1100, help='Maximum ' +
                        'number of active GCP VMs before skipping submission; ' +
                        'useful for throttling Cromwell server load')
    parser.add_argument('--contig-gate', type=int, default=len(hg38_primary_contigs), 
                        help='Maximum number of active contigs before skipping ' +
                        'submission; useful for throttling Cromwell server load')
    parser.add_argument('--submission-gate', type=float, default=0.5, help='Number ' +
                        'of minutes to pause after a new workflow submission ' +
                        '[default: 30 seconds]')
    parser.add_argument('--gcp-report-period', type=float, default=5, 
                        help='Number of minutes between GCP VM load reports ' +
                        '[default: report every 5 minutes]')
    parser.add_argument('--max-cycles', type=int, help='Maximum number of ' +
                        'workflow management cycles before exiting [default: ' +
                        'run until completion]')
    parser.add_argument('--max-runtime', type=int, help='Maximum number of ' +
                        'minutes to allow workflow management process to run ' +
                        'before exiting [default: run until completion]')
    parser.add_argument('--workflow-id-log-prefix', help='Optional prefix for ' +
                        'logger files for submitted workflow IDs')
    parser.add_argument('--dry-run', default=False, action='store_true',help='Do ' +
                        'not actually submit jobs or manipulate cloud objects.')
    parser.add_argument('--quiet', default=False, action='store_true',help='Do ' +
                        'not print logging, progress, and diagnostics to stdout')
    parser.add_argument('--verbose', default=False, action='store_true',
                        help='Print extra logging')
    args = parser.parse_args()

    # Infer method name if not provided
    wdl_name = path.splitext(path.basename(args.wdl))[0]
    if args.name is None:
        method_name = wdl_name
    else:
        method_name = args.name

    # Check to make sure --bucket is set
    bucket = args.bucket
    if args.bucket is None:
        bucket = getenv('WORKSPACE_BUCKET')
    if bucket is None:
        exit('Must provide a valid value for --bucket. ' + \
             'See --help for more information.')

    # Infer expected paths and other contig-specific constants
    input_json_fmt = sub('[/]+$', '', args.base_directory) + \
                     '/cromshell/inputs/' + \
                     sub('.template', '', 
                         path.splitext(path.basename(args.input_json_template))[0]) + \
                     '.{}.json'
    output_bucket_fmt = sub('[/]+$', '', args.staging_bucket) + '/{0}'
    output_json_fmt = output_bucket_fmt + '/' + method_name + '.{0}.outputs.json'
    if args.workflow_id_log_prefix is not None:
        prev_wid_fname_base = args.workflow_id_log_prefix + '.' + method_name
    else:
        prev_wid_fname_base = method_name
    prev_wids_fmt = sub('[/]+$', '', args.base_directory) + \
                    '/cromshell/job_ids/' + prev_wid_fname_base + \
                    '.{}.job_ids.list'

    # Read custom list of contigs, if optioned
    if args.contig_list is None:
        contigs = hg38_primary_contigs
    else:
        with open(args.contig_list) as clist:
            contigs = [x.rstrip() for x in clist.readlines()]
    contigs = g2cpy.chromsort(contigs)

    # Take the stricter requirement between --max-cycles and --max-runtime 
    # if either are specified
    max_cycles = []
    if args.max_cycles is not None:
        max_cycles.append(int(args.max_cycles))
    if args.max_runtime is not None:
        max_cycles.append(int(np.floor(args.max_runtime / args.outer_gate)))
    if len(max_cycles) > 0:
        max_cycles = np.nanmax([1, np.nanmin(max_cycles)])
    else:
        max_cycles = np.Inf

    # Report startup conditions
    if not args.quiet:
        with print_lock:
            startup_report({'Method name' : method_name,
                            'WDL' : args.wdl,
                            'Input .json template' : args.input_json_template,
                            'Cromwell options .json' : args.cromwell_options_json,
                            'WDL dependencies .zip' : args.dependencies_zip,
                            'Status tracker .tsv' : args.status_tsv,
                            'Workspace bucket' : bucket,
                            'Staging bucket' : args.staging_bucket,
                            'Local root' : args.base_directory,
                            'Behavior' : args.behavior,
                            'Dumpster' : args.dumpster,
                            'Automatic cleanup' : not args.no_cleanup,
                            'Hard reset' : args.hard_reset,
                            'Contigs' : ', '.join(contigs),
                            'Max attempts' : args.max_attempts,
                            'Outer gate' : args.outer_gate,
                            'Submission gate' : args.submission_gate,
                            'Active VM gate' : args.vm_gate,
                            'Active contig gate' : args.contig_gate,
                            'GCP load report period' : args.gcp_report_period,
                            'Max cycles' : max_cycles,
                            'Dry run' : args.dry_run,
                            'Quiet mode' : args.quiet})

    # If --status-tsv is provided and exists, load this file as a dict
    if args.status_tsv is not None \
    and path.exists(args.status_tsv):
        all_status = pd.read_csv(args.status_tsv, sep='\t', header=None).\
                        set_index(0)[1].to_dict()
        run_status = {k : v for k, v in all_status.items() if k in contigs}
    # Otherwise, create a dictionary for tracking status
    else:
        all_status = {k : 'unknown' for k in contigs}
        run_status = all_status.copy()

    # If --hard-reset is specified, unstage all contigs prior to entering monitor loop
    if args.hard_reset:

        # Check for staged outputs and update status if necessary
        # This is necessary if outputs are staged but no --status-tsv is provided
        for contig in contigs:
            if g2cpy.check_gcp_uris([output_json_fmt.format(contig)]):
                all_status[contig] = 'staged'
                run_status[contig] = 'staged'

        # Unstage staged outputs, if any exist
        if 'staged' in run_status.values():
            unstage_contigs = [k for k, v in run_status.items() if v == 'staged']
            for contig in unstage_contigs:
                if args.dry_run:
                    if not args.quiet:
                        msg = '[{}] Would un-stage outputs for {} if not for --dry-run'
                        with print_lock:
                            print(msg.format(clean_date(), contig))
                else:
                    msg = '[{}] Un-staging outputs for {}'
                    with print_lock:
                        print(msg.format(clean_date(), contig))
                    g2cpy.delete_uris(['/'.join([sub('/$', '', args.staging_bucket), 
                                                 contig, '**'])])
                    run_status[contig] = 'unstaged'
                    all_status.update(run_status)

    # Report beginning status
    if not args.quiet:
        with print_lock:
            report_status(run_status, title='Status at launch:')

    # Launch Cromwell load monitor
    if not args.quiet:
        t = threading.Thread(
            target=vm_load_daemon, 
            args=(args.gcp_report_period,), 
            daemon=True
            )
        t.start()

    # Loop infinitely while any status is not 'staged' and max_cycles has not been reached
    cycle_number = 0
    hard_reset_retries = {k : 0 for k in contigs}
    while cycle_number < max_cycles:
        cycle_number += 1
        if not args.quiet:
            msg = '[{}] Not all contigs yet staged. Entering another cycle of ' + \
                  'submission management routine.\n'
            with print_lock:
                print(msg.format(clean_date()))

        # Loop over all chromosomes
        for contig in contigs:

            # Get most recent status for this contig
            status = run_status.get(contig, 'unknown')

            # Format expected input + output paths
            input_json = input_json_fmt.format(contig)
            prev_wids = prev_wids_fmt.format(contig)
            output_bucket = output_bucket_fmt.format(contig)
            output_json_uri = output_json_fmt.format(contig)

            # Check possible statuses in a specific order to minimize necessary 
            # gsutil/cromshell queries
            while True:

                # If this is the first cycle and --hard-reset was specified,
                # always treat every submission as not_started
                if cycle_number == 1 and args.hard_reset:
                    status = 'not_started'
                    break

                # First, check if contig output has been already staged
                # If so, nothing more needs to be done
                if status == 'staged':
                    break
                if g2cpy.check_gcp_uris([output_json_uri]):
                    status = 'staged'
                    break

                # Second, check if contig has any previous workflow submissions
                # If not, contig has not yet been started and can be reported as such
                if status in 'not_started unknown'.split() \
                and not path.isfile(prev_wids):
                    status = 'not_started'
                    n_prior_subs = 0
                    break

                # If contig isn't either unstarted or fully staged, proceed with working 
                # through the possible intermediate states

                # Read ordered list of workflow IDs
                if path.exists(prev_wids):
                    with open(prev_wids) as fin:
                        wids = [line.rstrip() for line in fin.readlines()]
                        n_prior_subs = len(wids)
                else:
                    wids = []
                    n_prior_subs = 0

                # Check to make sure at least one workflow ID is reported
                if n_prior_subs == 0:
                    status = 'not_started'
                    break

                # Update status according to most recent workflow unless
                wid = wids[-1]
                status = g2cpy.check_workflow_status(wid, timeout=120)

                # If most recent workflow was successful, stage outputs and clear all 
                # files from Cromwell execution & output buckets
                if status == 'succeeded':
                    if args.dry_run:
                        if not args.quiet:
                            msg = '[{}] Found successful workflow ({}); would stage ' + \
                                  'outputs and clean garbage if not for --dry-run'
                            with print_lock:
                                print(msg.format(clean_date(), wid))
                    elif args.no_stage:
                        if not args.quiet:
                            msg = '[{}] Found successful workflow ({}); would stage ' + \
                                  'outputs and clean garbage if not for --no-stage'
                            with print_lock:
                                print(msg.format(clean_date(), wid))
                    else:
                        relocate_outputs(wid, output_bucket, wdl_name, output_json_uri)
                        status = 'staged'
                        if not args.no_cleanup:
                            g2cpy.collect_workflow_trash(wids, bucket, wdl_name,
                                                         args.dumpster)
                            subprocess.run('. code/refs/general_bash_utils.sh ' + \
                                           '&& cleanup_garbage', shell=True, check=False, 
                                           text=True, executable='/bin/bash')
                        if not args.quiet:
                            with print_lock:
                                print('')

                # Exit loop after processing most recent workflow ID
                break

            # For runs with --hard-reset, only count prior submissions
            # within the scope of this script
            if args.hard_reset:
                n_prior_subs = hard_reset_retries.get(contig, 0)

            # Submit new workflow if not already running or staged
            if status not in 'staged submitted running status_check_failed'.split():
                if n_prior_subs < args.max_attempts:

                    # Check to ensure that the cromwell server isn't already overloaded
                    n_active_vms = g2cpy.count_vms()
                    if n_active_vms > args.vm_gate:
                        if not args.quiet:
                            with print_lock:
                                msg1 = '[{}] Status of contig {}: {}'
                                print(msg1.format(clean_date(), contig, status))
                                msg2 = '[{}] Skipped submission of {} due to the ' + \
                                       'number of active VMs ({:,}) being greater ' + \
                                       'than the value of --vm-gate ({:,})'
                                print(msg2.format(clean_date(), contig, n_active_vms,
                                                  args.vm_gate))
                        status = 'vm_gate_skipped'

                    # Check to ensure that a permissible number of contigs are in flight
                    n_active_contigs = sum([s in running_status for s in run_status.values()])
                    if n_active_contigs > args.contig_gate:
                        if not args.quiet:
                            with print_lock:
                                msg1 = '[{}] Status of contig {}: {}'
                                print(msg1.format(clean_date(), contig, status))
                                msg2 = '[{}] Skipped submission of {} due to the ' + \
                                       'number of active contigs ({:,}) being greater ' + \
                                       'than the value of --contig-gate ({:,})'
                                print(msg2.format(clean_date(), contig, 
                                                  n_active_contigs,
                                                  args.contig_gate))
                        status = 'contig_gate_skipped'

                    # Otherwise, submit workflow as normal
                    else:
                        if args.behavior in 'indifferent strict'.split() \
                        or (args.behavior == 'patient' and status in strict_sub_status):
                            if args.behavior == 'strict':
                                aborted = g2cpy.abort_workflow(wid)
                            status = submit_workflow(contig, args.wdl, 
                                                     args.input_json_template, 
                                                     input_json, 
                                                     args.cromwell_options_json,
                                                     prev_wids, 
                                                     args.dependencies_zip, 
                                                     args.contig_variable_overrides, 
                                                     args.dry_run, args.quiet, 
                                                     args.verbose)
                            hard_reset_retries[contig] += 1
                            if args.submission_gate is not None:
                                if not args.quiet:
                                    msg = '[{}] Pausing {:,} minutes before ' + \
                                          'proceeding with next chromosome...\n'
                                    with print_lock:
                                        print(msg.format(clean_date(), 
                                                         args.submission_gate, 
                                                         args.submission_gate))
                                # Sleep for --submission-gate minutes
                                sleep(60 * args.submission_gate)

                else:
                    if not args.quiet:
                        msg = '[{}] Contig {} has reached the maximum number of ' + \
                              'attempts ({:,}) for {}; skipping.'
                        with print_lock:
                            print(msg.format(clean_date(), contig, 
                                             args.max_attempts, method_name))
                    status = 'exhausted'

            # Update & report status
            run_status[contig] = status
            all_status.update(run_status)
            if not args.quiet and status != 'vm_gate_skipped':
                msg = '[{}] Status of contig {}: {}'
                with print_lock:
                    print(msg.format(clean_date(), contig, status))
            if args.status_tsv is not None:
                with open(args.status_tsv, 'w') as fout:
                    for k, v in all_status.items():
                        fout.write('\t'.join([k, v]) + '\n')

        # Check at end of loop if all contigs are staged
        if len(set(run_status.values()).difference({'staged'})) == 0:
            break

        # Check if all contigs are exhausted; if so, no progress will be made
        if len(set(run_status.values()).difference({'staged', 'exhausted'})) == 0:
            msg = '[{}] All contigs are staged or exhausted. No more progress ' + \
                  'is possible with current settings. Exiting management routine.'
            with print_lock:
                print(msg.format(clean_date()))
            exit()

        # Wait before checking again
        else:
            if not args.quiet:
                with print_lock:
                    report_status(run_status)
                msg = '[{}] Finished checking status for all {:,} contigs. ' + \
                      'Waiting {:,} minutes before checking again...\n'
                with print_lock:
                    print(msg.format(clean_date(), len(contigs), args.outer_gate))

            # Sleep for --outer-gate minutes
            sleep(60 * args.outer_gate)

    # Report if time limit reached
    if cycle_number == max_cycles:
        msg = '[{}] Maximum number of cycles reached ({:,}). ' + \
              'Exiting management routine.'
        with print_lock:
            print(msg.format(clean_date(), max_cycles))
        exit(1)

    # Report completion once all contigs are staged
    msg = '[{}] All contig outputs have been staged! Exiting management routine.'
    with print_lock:
        print(msg.format(clean_date()))


if __name__ == '__main__':
    main()

