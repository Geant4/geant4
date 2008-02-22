#
# GEANT4 SBT Script to test G4Sphere
# I.Hrivnacova, IPN Orsay 23/01/2008 
#
test/maxPoints 1000
#
# --- sphere.a1.log
# No inner radius,  no phi segment, no theta segment
#
/solid/G4Sphere 0 1 0 360 0 180
/test/errorFileName  log/sphere.a1.log
/test/run
/voxel/errorFileName log/spherev.a1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/sphere.a2.log
/test/run
#
# --- sphere.b1.log
# With inner radius,  no phi segment, no theta segment
#
/solid/G4Sphere 0.5 1 0 360 0 180
/test/gridSizes 0 0 0 m
/test/errorFileName  log/sphere.b1.log
/test/run
/voxel/errorFileName log/spherev.b1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/sphere.b2.log
/test/run
#
# --- sphere.c1.log
# No inner radius,  with phi segment, no theta segment
#
/solid/G4Sphere 0 1 0 90 0 180
/test/gridSizes 0 0 0 m
/test/errorFileName  log/sphere.c1.log
/test/run
/voxel/errorFileName log/spherev.c1.log
/voxel/run
/test/gridSizes 0.02 0.02 0.02 m
/test/errorFileName  log/sphere.c2.log
/test/run
#
# --- sphere.d1.log
# With inner radius,  with phi segment, no theta segment
#
/solid/G4Sphere 0.5 1 0 90 0 180
/test/gridSizes 0 0 0 m
/test/errorFileName  log/sphere.d1.log
/test/run
/voxel/errorFileName log/spherev.d1.log
/voxel/run
/test/gridSizes 0.02 0.02 0.02 m
/test/errorFileName  log/sphere.d2.log
/test/run
#
# --- sphere.e1.log
# No inner radius,  no phi segment, with theta segment
#
/solid/G4Sphere 0 1 0 360 0 90
/test/gridSizes 0 0 0 m
/test/errorFileName  log/sphere.e1.log
/test/run
/voxel/errorFileName log/spherev.e1.log
/voxel/run
/test/gridSizes 0.02 0.02 0.02 m
/test/errorFileName  log/sphere.e2.log
/test/run
#
# --- sphere.f1.log
# With inner radius,  no phi segment, with theta segment
#
/solid/G4Sphere 0.5 1 0 360 0 90
/test/gridSizes 0 0 0 m
/test/errorFileName  log/sphere.f1.log
/test/run
/voxel/errorFileName log/spherev.f1.log
/voxel/run
/test/gridSizes 0.02 0.02 0.02 m
/test/errorFileName  log/sphere.f2.log
/test/run
#
# --- sphere.g1.log
# No inner radius,  with phi segment, with theta segment
#
/solid/G4Sphere 0 1 0 90 0 90
/test/gridSizes 0 0 0 m
/test/errorFileName  log/sphere.g1.log
/test/run
/voxel/errorFileName log/spherev.g1.log
/voxel/run
/test/gridSizes 0.02 0.02 0.02 m
/test/errorFileName  log/sphere.g2.log
/test/run
#
# --- sphere.h1.log
# With inner radius,  with phi segment, with theta segment
#
/solid/G4Sphere 0.5 1 0 90 0 90
/test/gridSizes 0 0 0 m
/test/errorFileName  log/sphere.h1.log
/test/run
/voxel/errorFileName log/spherev.h1.log
/voxel/run
/test/gridSizes 0.02 0.02 0.02 m
/test/errorFileName  log/sphere.h2.log
/test/run
#
# --- sphere.i1.log
# With inner radius,  with phi segment, with theta segment
# Very thin in radius
#
/solid/G4Sphere 0.00999 0.01001 10 260 30 70
/test/gridSizes 0.005 0.005 0.005 m
/test/errorFileName  log/sphere.i1.log
/test/run
/voxel/errorFileName log/spherev.i1.log
/voxel/run
exit
