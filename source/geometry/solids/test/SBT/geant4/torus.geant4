#
# GEANT4 SBT Script to test G4Torus
# I.Hrivnacova, IPN Orsay 23/01/2008 
#
/test/maxPoints 10000
#
# --- torus.a1.log
# No inner radius, no phi segment
#
/solid/G4Torus 0 0.4 1 0 360
/test/errorFileName log/torus.a1.log
/test/run
/voxel/errorFileName log/torusv.a1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/torus.a2.log
/test/run
#
# --- torus.b1.log
# With inner radius, no phi segment
#
/solid/G4Torus 0.2 0.4 1 0 360
/test/gridSizes 0 0 0 m
//test/errorFileName  log/torus.b1.log
/test/run
/voxel/errorFileName log/torusv.b1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/torus.b2.log
/test/run
#
# --- torus c1.log
# No inner radius, with phi segment
#
/solid/G4Torus 0 0.4 1 0 90
/test/gridSizes 0 0 0 m
/test/errorFileName  log/torus.c1.log
/test/run
/voxel/errorFileName log/torusv.c1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/torus.c2.log
/test/run
#
# --- torus d1.log
# With inner radius, with phi segment
#
/solid/G4Torus 0.2 0.4 1 0 90
/test/gridSizes 0 0 0 m
/test/errorFileName  log/torus.d1.log
/test/run
/voxel/errorFileName log/torusv.d1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/torus.d2.log
/test/run
#
# --- torus e1.log
# With very thin radius, with phi segment
#
/solid/G4Torus 0.399 0.4 1 0 90
/test/gridSizes 0 0 0 m
/test/errorFileName  log/torus.e1.log
/test/run
/voxel/errorFileName log/torusv.e1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/torus.e2.log
/test/run
#
exit
