#
# GEANT4 SBT Script to test G4Tubs
# DCW 19/3/99 First try
#
# Increment the number below to **really** waste CPU time
#
test/maxPoints 1000
#
# --- twistedTubs.a1.log
# Twisted tube with one phi segment
#
/test/gridSizes 0 0 0 m
/solid/G4TwistedTubs 30 0.8 1 -1 1 1 90
/test/errorFileName  log/twistedTubs.a1.log
/test/run
/voxel/errorFileName log/twistedTubsv.a1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/twistedTubs.a2.log
/test/run
#
exit
#
# --- twistedTubs.b1.log
# Twisted tube with more phi segments
# (visualization shows just one)
#
//solid/G4TwistedTubs 30 0.8 1 -1 1 9 90
/test/gridSizes 0 0 0 m
/test/errorFileName  log/twistedTubs.b1.log
/test/run
/voxel/errorFileName log/twistedTubsv.b1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/twistedTubs.b2.log
/test/run
#
exit
