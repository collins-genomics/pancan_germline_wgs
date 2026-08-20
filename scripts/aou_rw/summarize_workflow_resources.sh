#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Summarize resource usage for all tasks from single Cromwell workflow

usage() {
cat << EOF

USAGE: ./summarize_workflow_resources.sh workflow_id [execution_bucket]

       Collates data from montitoring logs from a workflow executed via Cromwell

       Options:
       - Can provide a file listing multiple workflow IDs instead of single workflow ID
       - Can provide execution bucket URI as a second argument
         (default bucket: \$WORKSPACE_BUCKET/workflows/cromwell-executions/)

EOF
}

# Parse arguments
wid_or_list=$1
EXEC_BUCKET=$2 # Optional prefix for Cromwell execution bucket
if [ $# -lt 1 ]; then
  echo -e "\nError: must supply workflow ID as first positional argument" >&2
  usage
  exit
fi

if [ -z $EXEC_BUCKET ]; then
  EXEC_BUCKET="$WORKSPACE_BUCKET/workflows/cromwell-executions"
fi

WRKDIR=`mktemp -d`

# Determine list of WIDs to process
if [ -s $wid_or_list ]; then
  # echo -e "Reading workflow IDs from $wid_or_list\n"
  cp $wid_or_list $WRKDIR/wids.list
else
  echo "$wid_or_list" > $WRKDIR/wids.list
fi

# Determine number of processors to use for parallelism
PROCS=$(( $(nproc) * 2 ))
if [ $PROCS -gt 64 ]; then
  PROCS=64
fi

# Loop over each workflow ID and process all logs in serial
echo -e "#task\tresource\tallocated\tpeak_used" > $WRKDIR/resources.tsv
while read wid; do

  # Make workflow-specific staging directory
  mkdir $WRKDIR/$wid
  SUBDIR=$WRKDIR/$wid

  # Gather list of monitoring.logs
  gsutil -m ls $EXEC_BUCKET/*/$wid/**monitoring.log 2>/dev/null \
  > $SUBDIR/log.uris.list 
  if [ $( cat $SUBDIR/log.uris.list | wc -l ) -lt 1 ]; then
    echo -e "\nWarning: no files named monitoring.log found in $EXEC_BUCKET/*/$wid/" >&2
    continue
  fi

  # Localize monitoring logs while mapping to directory structure of URIs
  n_logs_remote=$( cat $SUBDIR/log.uris.list | wc -l )
  # echo -e "Identified $n_logs_remote logs for workflow $wid. Now localizing.\n"
  while read uri; do
    # Recursively create all subdirectories
    local_log_dir=$SUBDIR/$( dirname $( echo $uri | sed "s/$wid\//\t/g" | awk '{ print $2 }' ) )
    mkdir -p $local_log_dir
    
    # Print source and destination to stream to gsutil cp -I
    echo -e "$uri\t$local_log_dir/monitoring.log"
  done < $SUBDIR/log.uris.list \
  > $SUBDIR/log.loc.map.tsv
  cat $SUBDIR/log.loc.map.tsv \
  | xargs -n 2 -P $PROCS sh -c 'gsutil -m cp "$1" "$2"' _ 2>/dev/null
  n_logs_local=$( find $SUBDIR/ -name "*monitoring.log" | wc -l )
  if [ $n_logs_remote -ne $n_logs_local ]; then
    echo -e "\nWarning: only successfully localized $n_logs_local of $n_logs_remote logs for $wid\n" >&2
  fi

  # Extract resource information from each monitoring log
  SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
  while read uri loc_log; do
    $SCRIPT_DIR/parse_cromwell_monitoring_log.py $loc_log
  done < $SUBDIR/log.loc.map.tsv \
  >> $WRKDIR/resources.tsv

  # Clean up intermediate garbage
  rm -rf $SUBDIR

done < $WRKDIR/wids.list

# Summarize resources across all workflows by task & print to stdout
$SCRIPT_DIR/summarize_workflow_resource_distribs.R $WRKDIR/resources.tsv

# Clean up
rm -rf $WRKDIR
