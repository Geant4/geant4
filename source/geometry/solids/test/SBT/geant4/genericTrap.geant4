#
# GEANT4 SBT Script to test G4GenericTrap
# derived from G4ExtrudedSolid done by I.Hrivnakova
# 19/11/2009 T.Nikitina

#
/test/maxPoints 1000
#
# --- genericTrap.a1.log
# Generic Trap with 8 vertices no twist(box like)
#

/solid/G4GenericTrap 1 (-3, -3, 3, 3, -3, -3, 3, 3) (-3, 3, 3, -3, -3, 3, 3, -3)
/test/gridSizes 0.1 0.1 0.1 m
/test/errorFileName  log/genericTrap.a1.log
/test/run
/voxel/errorFileName log/genericTrapv.a1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/genericTrap.a2.log
/test/run
#
# --- genericTrap.b1.log
# Generic Trap with 8 vertices no twist(Trd like)
#

/solid/G4GenericTrap 1 (-3, -3, 3, 3, -3, -3, 3, 3) (-2, 2, 2, -2, -2, 2, 2, -2)
/test/gridSizes 0.1 0.1 0.1 m
/test/errorFileName  log/genericTrap.b1.log
/test/run
/voxel/errorFileName log/genericTrapv.b1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/genericTrap.b2.log
/test/run
#
# --- genericTrap.c1.log
# Generic Trap with 5 dif. vertices no twist(Tet like)
#

/solid/G4GenericTrap 1 (-3, -3, 3, 3, -3, -3, 3, 3) (-2, -2, -2, -2, -2, -2, -2, -2)
/test/gridSizes 0.1 0.1 0.1 m
/test/errorFileName  log/genericTrap.c1.log
/test/run
/voxel/errorFileName log/genericTrapv.c1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/genericTrap.c2.log
/test/run
#
#
# --- genericTrap.d1.log
# Generic Trap with 8 vertices with twist
#

/solid/G4GenericTrap 1 (-3, -3, 3, 3, -0.5, -2, 2, 2) (-3, 3, 3, -3, -2, 2, 2, -2)
/test/gridSizes 0.1 0.1 0.1 m
/test/errorFileName  log/genericTrap.d1.log
/test/run
/voxel/errorFileName log/genericTrapv.d1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/genericTrap.d2.log
/test/run
#
exit
