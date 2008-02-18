#
# GEANT4 SBT Script to test G4Box
# I.Hrivnacova, IPN Orsay 29/01/2008 
# 
#
/test/maxPoints 10000
#
# --- twistedBox.a1.log
# Box with a twisted angle
#
/solid/G4TwistedBox 30 1 1 1
/test/errorFileName log/twistedBox.a1.log
/test/run
/voxel/errorFileName log/twistedBoxv.a1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/twistedBox.a2.log
/test/run
#
# As the test is very time consuming
# we skip the other variants for the moment
exit
#
# --- twistedBox.b1.log
# Box with a very small twisted angle
#
/solid/G4TwistedBox 1 1 1 1
/test/gridSizes 0 0 0 m
/test/errorFileName log/twistedBox.b1.log
/test/run
/voxel/errorFileName log/twistedBoxv.b1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/twistedBox.b2.log
/test/run
#
# --- twistedBox.c1.log
# Very thin box with a twisted angle
#
/solid/G4TwistedBox 30 0.0001 1 2
/test/gridSizes 0 0 0 m
/test/errorFileName log/twistedBox.c1.log
/test/run
/voxel/errorFileName log/twistedBoxv.c1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/twistedBox.c2.log
/test/run
#
# --- twistedBox.d1.log
# Very thin box with a very small twisted angle
#
/solid/G4TwistedBox 1 0.0001 1 2
/test/gridSizes 0 0 0 m
/test/errorFileName log/twistedBox.d1.log
/test/run
/voxel/errorFileName log/twistedBoxv.d1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/twistedBox.d2.log
/test/run
#
exit
