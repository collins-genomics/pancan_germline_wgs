# The Germline Genomics of Cancer (G2C)
## Variant calling from germline whole-genome sequencing (WGS) across cancer types

Copyright (c) 2023-Present, [Ryan L. Collins](mailto:Ryan_Collins@dfci.harvard.edu), the [Collins Genomics Lab](https://labs.dana-farber.org/collins-genomics/), and the [Van Allen Lab](https://labs.dana-farber.org/vanallenlab/) at the Dana-Farber Cancer Institute.  

Distributed under terms of the [GNU GPL v2.0 License](/LICENSE) (see `LICENSE`).  

#### _Note: this repository is under active development. More documentation will be added as the project evolves._

---  

## Synopsis    

This repository contains the working code and scripts used to detect, genotype, filter, and annotate germline variants from germline WGS across cancer types  

---  

## Table of Contents  

| Directory | Description |  
| :--- | :--- |  
| [`docker/`](https://github.com/talkowski-lab/dsmap/tree/main/docker) | Instructions for building project-related Docker images |   
| [`refs/`](https://github.com/talkowski-lab/dsmap/tree/main/refs) | Reference .jsons and dotfiles |   
| [`scripts/`](https://github.com/talkowski-lab/dsmap/tree/main/scripts) | Stand-alone scripts called by various workflows |   
| [`shell/`](https://github.com/talkowski-lab/dsmap/tree/main/shell) | Shell snippets for running specific processes or analyses |  
| [`src/`](https://github.com/talkowski-lab/dsmap/tree/main/src) | Source code for the `g2cpy` and `G2CR` companion libraries |  
| [`wdl/`](https://github.com/talkowski-lab/dsmap/tree/main/wdl) | Stand-alone WDL workflows |   

---  

## Data Access  

All site-frequency information and derived summary statistics will be made publicly available from this project to the extent permissible by sample consent and data use agreements. At this time, no public data is yet available, and we do not anticipate any public data becoming available prior to 2026.  

All sample-level data is stored in a secure Google Cloud bucket. Note that permissions must be granted for bucket access, and permissions cannot be granted to any users outside of those explicitly covered by the numerous data use agreements between DFCI/Van Allen lab and other organization. Contact [Ryan L. Collins](mailto:Ryan_Collins@dfci.harvard.edu) or [Erin Shannon](ErinE_Shannon@dfci.harvard.edu) if you believe you qualify for access.  
