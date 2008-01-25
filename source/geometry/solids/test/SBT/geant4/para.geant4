#
# GEANT4 SBT Script to test G4Para
# I.Hrivnacova, IPN Orsay 23/01/2008 
#
/test/maxPoints 10000
#
# --- para.a1.log
# No alpha, no theta, no phi
#
/solid/G4Para 1 1 1 0 0 0
/test/errorFileName log/para.a1.log
/test/run
/voxel/errorFileName log/parav.a1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/para.a2.log
/test/run
#
# --- para.b1.log
# With alpha, no theta, no phi
#
/solid/G4Para  1 1 1 30 0 0
/test/gridSizes 0 0 0 m
/test/errorFileName  log/para.b1.log
/test/run
/voxel/errorFileName log/parav.b1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/para.b2.log
/test/run
#
# --- para c1.log
# With alpha, with theta, no phi
#
/solid/G4Para  1 1 1 30 30 0
/test/gridSizes 0 0 0 m
/test/errorFileName  log/para.c1.log
/test/run
/voxel/errorFileName log/parav.c1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/para.c2.log
/test/run
#
# --- para d1.log
# With alpha, with theta, with phi
#
/solid/G4Para  1 1 1 30 30 30
/test/gridSizes 0 0 0 m
/test/errorFileName  log/para.d1.log
/test/run
/voxel/errorFileName log/parav.d1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/para.d2.log
/test/run
#
# --- para e1.log
# With alpha, with theta, with phi
# Make non equal sides, one very thin
#
/solid/G4Para  0.001 1 2 30 30 30
/test/gridSizes 0 0 0 m
/test/errorFileName  log/para.e1.log
/test/run
/voxel/errorFileName log/parav.e1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/para.e2.log
/test/run
#
exit
