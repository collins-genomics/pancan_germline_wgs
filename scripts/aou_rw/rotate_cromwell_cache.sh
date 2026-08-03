#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Helper script to rotate Cromwell cache

DBDIR=/home/jupyter/.cromwell/db

if [[ -d "$DBDIR" ]]; then

    size_gb=$( du -sBG "$DBDIR" | cut -f1 | tr -d 'G' )

    if [ $# -gt 0 ] && [ $1 == "delete" ]; then

        echo "================================"
        echo "Clearing Cromwell cache (${size_gb} GB)"
        echo "================================"

        rm -rf "$DBDIR"

    else

        echo "================================"
        echo "Rotating Cromwell cache (${size_gb} GB)"
        echo "================================"

        mv -f "$DBDIR" "${DBDIR}.previous"

    fi

fi

mkdir -p "$DBDIR"