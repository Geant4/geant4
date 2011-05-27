#
# GEANT4 SBT Script to test G4GenericTrap
# derived from G4ExtrudedSolid done by I.Hrivnakova
# 19/11/2009 T.Nikitina

#
/test/maxPoints 1000
#
# --- extrudedSolid.a1.log
# Generic Trap with 8 vertices
#
/solid/G4Paraboloid 1 1 2
/test/gridSizes 0.1 0.1 0.1 m
/test/errorFileName  log/paraboloid.a1.log
/test/run
/voxel/errorFileName log/paraboloidv.a1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/paraboloid.a2.log
/test/run
#
exit
