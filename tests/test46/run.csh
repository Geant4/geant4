#/bin/csh

mkdir -p $REFERENCE
cd $REFERENCE
mkdir -p $PHYSLIST
cd $PHYSLIST

set part = "pi-"

source run_part.csh ${part}

echo "Done!"
#
