#/bin/csh
#
# Run single configuration
#

mkdir -p $REFERENCE
cd $REFERENCE
mkdir -p $PHYSLIST
cd $PHYSLIST

source $G4INSTALL/tests/test46/run_single.csh $1

cd ../
echo $REFERENCE/$PHYSLIST "Done!"
#
