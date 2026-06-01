#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell subroutine to install common dependencies not preinstalled on AoU nodes

# Note that this code is designed to be run inside the AoU Researcher Workbench

# Any positional arguments given to this script will be used to specify which
# subsets of libraries (i.e., which languages) should be installed
if [ $# -eq 0 ]; then
  echo -e "Error: at least one positional argument must be supplied to install_packages.sh"
  echo -e "Current options include: R, python"
  echo -e "Exiting"
  exit 1
fi

# Check to ensure this is executed from the correct directory
for dir in code code/src; do
  if ! [ -e $dir ]; then
    echo -e "Error: could not locate $dir directory."
    echo -e "Are you sure you are calling this script from the same directory where G2C code has been staged?"
    echo -e "Exiting"
    exit 1
  fi
done

# Install all packages as optioned
for lang in "$@"; do

  case $lang in

    # Install R libraries
    R)

      # Install all R libraries distributed via CRAN
      for lib in argparse optparse beeswarm bedr caret EQL vioplot DescTools colourvalues viridis; do
        Rscript -e "if(require('$lib') == FALSE){install.packages('$lib', repos='https://cloud.r-project.org')}"
      done

      # Install RLCtools
      export RLCtools_version=0.1
      Rscript -e "if(require('RLCtools') == TRUE){remove.packages('RLCtools')}; install.packages('code/src/RLCtools_$RLCtools_version.tar.gz', repos=NULL, type='source')"

      # Install G2C companion library
      export G2CR_version=0.3.0
      Rscript -e "if(require('G2CR') == TRUE){remove.packages('G2CR')}; install.packages('code/src/G2CR_$G2CR_version.tar.gz', repos=NULL, type='source')"
      ;;

    # Install python packages
    python)
      
      # Ensure pip version 25.0 is used for compatability with editable install for svtk
      pip install pip==25.0

      # Install various public python packages
      for pkg in Cython; do
        pip install $pkg
      done

      # Install svtk & G2C companion package from source
      pip install -e code/src/g2cpy
      pip install -e code/src/svtk
      ;;

  esac
done
