#
# GEANT4 SBT Script to test G4Tet
# I.Hrivnacova, IPN Orsay 28/01/2008 

#
test/maxPoints 1000
#
# --- tet.a1.log
# Regular
#
/solid/G4Tet (0.0,0.0,1.73205080756887719) (0.0,1.63299316185545207,-0.577350269189625842) (-1.41421356237309515,-0.816496580927726034,-0.577350269189625842) (1.41421356237309515,-0.816496580927726034,-0.577350269189625842) 
/test/gridSizes 0 0 0 m
/test/errorFileName  log/tet.a1.log
/test/run
/voxel/errorFileName log/tetv.a1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/tet.a2.log
/test/run
#
# --- tet.b1.log
# Assymetric
#
/solid/G4Tet (0.0,0.0,1.0) (-1.0,-1.0,-1.0) (+1.0,-1.0,-1.0) (0.0,1.0,-1.0)
/test/gridSizes 0 0 0 m
/test/errorFileName  log/tet.b1.log
/test/run
/voxel/errorFileName log/tetv.b1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/tet.b2.log
/test/run
#
# --- tet.c1.log
# Assymetric, more extreme
#
/solid/G4Tet (0.0,0.0,1.0) (-0.2,-1.0,-1.0) (0.2,-1.0,-1.0) (0.0,1.0,-1.0)
/test/gridSizes 0 0 0 m
/test/errorFileName  log/tet.c1.log
/test/run
/voxel/errorFileName log/tetv.c1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/tet.c2.log
/test/run
#
# --- tet.d1.log
# Assymetric, very extreme
#
/solid/G4Tet (0.0,0.0,1.0) (-0.001,-1.0,-1.0) (+0.001,-1.0,-1.0) (0.0,1.0,-1.0)
/test/gridSizes 0 0 0 m
/test/errorFileName  log/tet.d1.log
/test/run
/voxel/errorFileName log/tetv.d1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/tet.d2.log
/test/run
#
exit
