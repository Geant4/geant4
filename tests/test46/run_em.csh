#/bin/csh
#
# Run single configuration
#
mkdir -p $PHYSLIST
cd $PHYSLIST

source $G4INSTALL/tests/test46/run_single.csh $1

cd ../
#
