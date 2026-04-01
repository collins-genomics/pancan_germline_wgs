#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Setup script for G2C python helper library
"""

from setuptools import find_packages, setup

setup(
    name='g2cpy',
    version='0.1.0',
    packages=['g2cpy'],
    author='Ryan L. Collins',
    author_email='Ryan_Collins@dfci.harvard.edu',
    url='https://github.com/collins-genomics/pancan_germline_wgs',
    description='Dana-Farber Germline Genomics of Cancer Database: python helper functions',
    long_description=open('README.md').read(),
    license='GNU GPL v2.0',
    install_requires=[],
)