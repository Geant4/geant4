#
# GEANT4 SBT Script to test G4Cons
# I.Hrivnacova, IPN Orsay 30/01/2008 
#
# Increment the number below when errors become a bit more rare
#
test/maxPoints 10000
#
# --- twistedTrap.a1.log
# Here is a trap,  used and getting problems in CMS    May 28, 1999
# + now plus twisted angle
#
/solid/G4TwistedTrap2   30 1.268 0 0 0.295 1.7122 1.87029  0.295  1.7122 1.87029 0
/test/gridSizes 0 0 0 m
/test/errorFileName  log/twistedTrap.a1.log
/test/run
/voxel/errorFileName log/twistedTrapv.a1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/twistedTrap.a2.log
/test/run
#
exit

#
# --- twistedTrap.b1.log
# Very easy case
#
/solid/G4TwistedTrap 30 1 1 1 1 
/test/gridSizes 0 0 0 m
/test/widths 1 1 1 m
/test/errorFileName log/twistedTrap.b1.log
/test/run
/voxel/errorFileName log/twistedTrapv.b1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/twistedTrap.b2.log
/test/run
#
# --- twistedTrap.c1.log
# Adjust just x
#
/solid/G4TwistedTrap 30 0.5 1.5 1 1
/test/gridSizes 0 0 0 m
/test/errorFileName  log/twistedTrap.c1.log
/test/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/twistedTrap.c2.log
/test/run
#
# --- twistedTrap.d1.log
# Extreme case
#
/solid/G4TwistedTrap 30 0.000001 1 0.00002 1
/test/widths 1 0.00002 1 m
/test/gridSizes 0 0 0 m
/test/errorFileName  log/twistedTrap.d1.log
/test/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/twistedTrap.d2.log
/test/run
#
exit
