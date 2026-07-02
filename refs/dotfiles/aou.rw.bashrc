#!/usr/bin/env bash

# .bash_profile for AoU RW v2.0 (Verily Pre)
# Updated June 2026
# Ryan Collins <Ryan_Collins@dfci.harvard.edu>

# Simple ls with options
alias l='ls -lhtr'

# Default to zless instead of less
alias 'less=zless'

# Colorize grep output
alias 'grep=grep --color=auto'
alias 'egrep=egrep --color=auto'
alias 'fgrep=fgrep --color=auto'

# Set timezone to Eastern US time
export TZ="America/New_York"
