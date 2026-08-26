#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Cromwell server launch script for AoU RW v2.0 (Verily Pre)

# Note that this must be run in a separate terminal

# Recommended VM configuration: 16 vCPU x 64GB RAM

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://rw-migration-aou-rw-84a0039b

# Source .bashrc and bash utility functions
. code/refs/dotfiles/aou.rw.bashrc
. code/refs/general_bash_utils.sh

# Set terminal timezone to Eastern US time
export TZDIR=$(python -c "import tzdata, os; print(os.path.join(os.path.dirname(tzdata.__file__), 'zoneinfo'))")
export TZ=America/New_York
date +"%Y-%m-%d %H:%M:%S %Z %z"

# Infer workspace number and save as environment variable
export WN=$( get_workspace_number )

# Reassign $WORKSPACE_BUCKET to match AoU v1-migrated URIs
case "$WN" in
  "w1")
    export WORKSPACE_BUCKET="gs://rw-migration-aou-rw-84a0039b"
    ;;
  "w2")
    export WORKSPACE_BUCKET="gs://rw-migration-aou-rw-78e2871d"
    ;;
  "w3")
    export WORKSPACE_BUCKET="gs://rw-migration-aou-rw-3c78b3b7"
    ;;
  "w4")
    export WORKSPACE_BUCKET="gs://rw-migration-aou-rw-efb2fd38"
    ;;
  "w5")
    export WORKSPACE_BUCKET="gs://rw-migration-aou-rw-484d2a66"
    ;;
  "dev")
    export WORKSPACE_BUCKET="gs://rw-migration-aou-rw-e34d8d8a"
    ;;
  *)
    echo "UNKNOWN WORKSPACE NUMBER"
    ;;
esac

# Create cromwell config file
if ! [ -e /home/jupyter/.cromwell ]; then
  mkdir /home/jupyter/.cromwell
fi
if [ -e /home/jupyter/.cromwell/cromwell.conf ]; then
  rm /home/jupyter/.cromwell/cromwell.conf
  rm /home/jupyter/.cromwell/cromwell.override.conf
fi
wb cromwell generate-config \
  --google-bucket-name=$WORKSPACE_BUCKET \
  --dir=/home/jupyter/.cromwell
cat << EOF > /home/jupyter/.cromwell/cromwell.override.conf
include "cromwell.conf"

backend.providers.GCPBATCH.config {
  concurrent-job-limit = 3000
}

call-caching {
  enabled = true
}

system {
  job-rate-control {
    jobs = 10
    per = 10 seconds
  }
  workflow-heartbeats {
    ttl = 20 minutes
    write-failure-shutdown-duration = 15 minutes
    write-batch-size = 250
  }
}

engine {
  filesystems {
    gcs {
      auth = "application_default"
    }
  }
}

database {
  profile = "slick.jdbc.HsqldbProfile$"

  db {
    driver = "org.hsqldb.jdbcDriver"

    url = """
      jdbc:hsqldb:file:/home/jupyter/.cromwell/db/cromwell;
      shutdown=false;
      hsqldb.default_table_type=cached;
      hsqldb.tx=mvcc;
      hsqldb.result_max_memory_rows=10000;
      hsqldb.large_data=true;
      hsqldb.script_format=3
    """

    connectionTimeout = 1200000
    numThreads = 8
  }

  insert-batch-size = 2000
}
EOF

# Ensure cromwell local database exists
DBDIR=/home/jupyter/.cromwell/db
mkdir -p $DBDIR
if [[ -d "$DBDIR" ]]; then
  size_gb=$( du -sBG "$DBDIR" | cut -f1 | tr -d 'G' )
  echo "Cromwell DB size: $size_gb GB"
fi

# Create cromshell config file
if [ ! -e /home/jupyter/.cromshell ]; then
  mkdir /home/jupyter/.cromshell
fi
cat << EOF > /home/jupyter/.cromshell/cromshell_config.json
{
  "cromwell_server": "http://localhost:8000",
  "requests_timeout": 5
}
EOF

# Launch cromwell in server mode
echo -e "\nNow launching cromwell server using $CROMWELL_JAR"
java -jar $CROMWELL_JAR --version
java \
  -Xms12G \
  -Xmx48G \
  -XX:+UseG1GC \
  -Dconfig.file=/home/jupyter/.cromwell/cromwell.override.conf \
  -Xlog:gc*:file=/home/jupyter/.cromwell/gc.log:time \
  -XX:+HeapDumpOnOutOfMemoryError \
  -XX:HeapDumpPath=/home/jupyter/.cromwell/ \
  -jar $CROMWELL_JAR \
  server
