#!/usr/bin/python

#-----------------------------------------------------------------
# This Python script is the main one. It has no input parameters.
# Look at the string "***LOOKHERE***" to see where you may want
# to change something, i.e. to specify which cases to run.
# This program assumes the existence of the two executables:
#     bin/Linux-g++/mainStatAccepTest-$REF1-$PHYSICS
#     bin/Linux-g++/mainStatAccepTest-$REF2-$PHYSICS
# for the Geant4 reference and Physics List considered.
# This Python script invokes the shell script  simuDriver.sh .
#-----------------------------------------------------------------

import os
import sys
import string

print '========== START mainScript.py =========='

#***LOOKHERE***
REF1        = "4.6.2.ref01"
REF2        = "4.6.2.ref03"
PHYSICS     = "LHEP"
EVENTS      = "100"
#***endLOOKHERE***

executable = "bin/Linux-g++/mainStatAccepTest-" + REF1 + "-" + PHYSICS
if ( not os.path.exists( executable ) ) :
    print '***ERROR*** : executable ', executable, '  NOT found!'
    sys.exit(0)
    
executable = "bin/Linux-g++/mainStatAccepTest-" + REF2 + "-" + PHYSICS
if ( not os.path.exists( executable ) ) :
    print '***ERROR*** : executable ', executable, '  NOT found!'
    sys.exit(0)

# ---------------------------------------------

#***LOOKHERE***
CALORIMETER = "FeSci"
PARTICLE    = "p"
ENERGY      = "10GeV"
#***endLOOKHERE***

os.system( "simuDriver.sh " + REF1 + " " + REF2 + " " + PHYSICS + " "
           + CALORIMETER + " " + PARTICLE + " " + ENERGY + " " + EVENTS )

#***LOOKHERE***
CALORIMETER = "CuLAr"
PARTICLE    = "pi+"
ENERGY      = "20GeV"
#***endLOOKHERE***

os.system( "simuDriver.sh " + REF1 + " " + REF2 + " " + PHYSICS + " "
           + CALORIMETER + " " + PARTICLE + " " + ENERGY + " " + EVENTS )

# ---------------------------------------------

print '========== END mainScript.py =========='
