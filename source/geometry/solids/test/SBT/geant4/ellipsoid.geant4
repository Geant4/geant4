#
# GEANT4 SBT Script to test G4Ellipsoid
# I.Hrivnacova, IPN Orsay 28/01/2008 
#
test/maxPoints 1000
#
# --- ellipsoid.a1.log
# Spheric, no z-cuts
#
/solid/G4Ellipsoid 1 1 1 0 0
/test/errorFileName  log/ellipsoid.a1.log
/test/run
/voxel/errorFileName log/ellipsoidv.a1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/ellipsoid.a2.log
/test/run
#
# --- ellipsoid.b1.log
# Elliptic, no z-cuts 
#
/solid/G4Ellipsoid 0.5 0.8 1 0 0
/test/gridSizes 0 0 0 m
/test/errorFileName  log/ellipsoid.b1.log
/test/run
/voxel/errorFileName log/ellipsoidv.b1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/ellipsoid.b2.log
/test/run
#
# --- ellipsoid.c1.log
# Elliptic,  with bottom z cut
#
/solid/G4Ellipsoid 0.5 0.8 1 -0.4 10
//test/gridSizes 0 0 0 m
/test/errorFileName  log/ellipsoid.c1.log
/test/run
/voxel/errorFileName log/ellipsoidv.c1.log
/voxel/run
/test/gridSizes 0.02 0.02 0.02 m
/test/errorFileName  log/ellipsoid.c2.log
/test/run
#
# --- ellipsoid.d1.log
# Elliptic,  with top z cut
#
/solid/G4Ellipsoid 0.5 0.8 1 -10 0.8
/test/gridSizes 0 0 0 m
/test/errorFileName  log/ellipsoid.d1.log
/test/run
/voxel/errorFileName log/ellipsoidv.d1.log
/voxel/run
/test/gridSizes 0.02 0.02 0.02 m
/test/errorFileName  log/ellipsoid.d2.log
/test/run
#
# --- ellipsoid.e1.log
# Elliptic,  with both bottom and top z cuts
#
/solid/G4Ellipsoid 0.5 0.8 1 -0.4 0.8
/test/gridSizes 0 0 0 m
/test/errorFileName  log/ellipsoid.e1.log
/test/run
/voxel/errorFileName log/ellipsoidv.e1.log
/voxel/run
/test/gridSizes 0.02 0.02 0.02 m
/test/errorFileName  log/ellipsoid.e2.log
/test/run
exit
