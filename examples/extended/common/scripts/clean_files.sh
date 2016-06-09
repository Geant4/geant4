#!/bin/sh

# Script to remove classes copied from common repository in examples.
# All files with ExG4 prefix will be removed.
#
# Usage:
# clean_files.sh [--force]

FROM=`pwd`

FORCE=""
if [ "$1" = "--force" ]; then
  FORCE="-f"
fi

rm $FORCE include/ExG4*.hh
rm $FORCE src/ExG4*.cc




