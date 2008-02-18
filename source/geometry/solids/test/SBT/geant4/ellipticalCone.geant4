#
# GEANT4 SBT Script to test G4EllipticalCone
# I.Hrivnacova, IPN Orsay 28/01/2008 

#
test/maxPoints 1000
#
# --- ellipticalCone.a1.log
# Symmetric, cutted
#
/solid/G4EllipticalCone 0.3 0.3 0.75 0.25
/test/errorFileName  log/ellipticalCone.a1.log
/test/run
/voxel/errorFileName log/ellipticalConev.a1.log
/voxel/run
#
# !!! test with grid causes exception:
# ZMxpvInfiniteVector thrown:
# ZMxpvInfiniteVector: Attempt to do vector /= 0 -- division by zero would produce infinite or NAN components
# at line 345 in file ThreeVector.cc
# terminate called after throwing an instance of 'ZMxpvInfiniteVector'
# what():  ZMxpvInfiniteVector: Attempt to do vector /= 0 -- division by zero would produce infinite or NAN components
#/test/gridSizes 0.2 0.2 0.2 m
#/test/errorFileName  log/ellipticalCone.a2.log
#/test/run
#
# --- ellipticalCone.b1.log
# Symmetric, with a vertex
#
/solid/G4EllipticalCone 0.3 0.3 0.75 0.75
/test/errorFileName  log/ellipticalCone.b1.log
/test/run
/voxel/errorFileName log/ellipticalConev.b1.log
/voxel/run
#/test/gridSizes 0.2 0.2 0.2 m
#/test/errorFileName  log/ellipticalCone.b2.log
#/test/run
#
# --- ellipticalCone.c1.log
# Elliptic, cutted
#
/solid/G4EllipticalCone 0.3 0.6 0.75 0.25
/test/errorFileName  log/ellipticalCone.c1.log
/test/run
/voxel/errorFileName log/ellipticalConev.c1.log
/voxel/run
#/test/gridSizes 0.2 0.2 0.2 m
#/test/errorFileName  log/ellipticalCone.c2.log
#/test/run
#
# --- ellipticalCone.d1.log
# Elliptic, with a vertex
#
/solid/G4EllipticalCone 0.3 0.6 0.75 0.75
/test/errorFileName  log/ellipticalCone.d1.log
/test/run
/voxel/errorFileName log/ellipticalConev.d1.log
/voxel/run
#/test/gridSizes 0.2 0.2 0.2 m
#/test/errorFileName  log/ellipticalCone.d2.log
#/test/run
exit
