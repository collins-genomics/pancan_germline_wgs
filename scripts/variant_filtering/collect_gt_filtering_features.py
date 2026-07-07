#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Gather genotype-level annotations for GT filtering model
Input VCF must strictly contain only a single variant type (snv, indel, or SVs)
"""


import argparse
import gzip
import g2cpy
import pysam
from sys import stdin, stdout


