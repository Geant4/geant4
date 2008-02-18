#
# GEANT4 SBT Script to test G4Cons
# I.Hrivnacova, IPN Orsay 31/01/2008 
#
# Increment the number below when errors become a bit more rare
#
test/maxPoints 10000
#
# --- twistedTrd.c1.log
# Adjust x and y
#
/solid/G4TwistedTrd 0.5 1.5  0.25 1 1 30 
/test/gridSizes 0 0 0 m
/test/errorFileName  log/twistedTrd.c1.log
/test/run
/voxel/errorFileName log/twistedTrdv.c1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/twistedTrd.c2.log
/test/run
#
exit
#
# --- twistedTrd.a1.log
# Very easy case
#
/solid/G4TwistedTrd 1 1 1 1 1 30 
/test/gridSizes 0 0 0 m
/test/widths 1 1 1 m
/test/errorFileName log/twistedTrd.a1.log
/test/run
/voxel/errorFileName log/twistedTrdv.a1.log
/voxel/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName log/twistedTrd.a2.log
/test/run
#
# --- twistedTrd.b1.log
# Adjust just x
#
/solid/G4TwistedTrd 0.5 1.5 1 1 1 30 
/test/gridSizes 0 0 0 m
/test/errorFileName  log/twistedTrd.b1.log
/test/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/twistedTrd.b2.log
/test/run
#
# --- twistedTrd.d1.log
# Extreme case
#
/solid/G4TwistedTrd 0.000001 1 0.00001 0.00002 1 30
/test/widths 1 0.00002 1 m
/test/gridSizes 0 0 0 m
/test/errorFileName  log/twistedTrd.d1.log
/test/run
/test/gridSizes 0.2 0.2 0.2 m
/test/errorFileName  log/twistedTrd.d2.log
/test/run
#
exit
