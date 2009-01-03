#/bin/csh

mkdir -p $REFERENCE
cd $REFERENCE
mkdir -p $PHYSLIST
cd $PHYSLIST

source $G4INSTALL/tests/test46/run_part.csh pi-

cd ../../
echo $REFERENCE/$PHYSLIST "Done!"
#
