#!/usr/bin/python

#-----------------------------------------------------------------
# This Python script is the main one. It has no input parameters.
# Look at the string "***LOOKHERE***" to see where you may want
# to change something, i.e. to specify which cases to run.
# This Python script invokes the shell script  simuDriver.sh .
#-----------------------------------------------------------------

import os
import sys
import string

print '========== START mainScript.py =========='

#***LOOKHERE*** You have to specify:
#               -  the two references to be compared:
#                  eg:   "6.2.p02"  and  "7.0.cand01"
#               -  whether the Geant4 simulation should be run for
#                  the two reference:
#                  "Yes" when you want that the simulation is run;
#                  "No"  when you do not want that the simulation is run.
#               -  whether you want to run the statistical tests:
#                  "Yes" when you want to run the statistical tests;
#                  "No"  when you do not want to run the statistical tests.
#               -  the Physics Lists:
#                  "LHEP", "QGSP", "QGSC", "QGSP_BERT", "QGSP_BIC", "QGSP_GN".
#               -  the Number of Events:
#                  e.g. : "5000", "5k", etc.
#               -  the Bfield:
#                  e.g. : "0", "4tesla", etc.
#
REF1        = "9.2.p02"
SIM_REF1    = "Yes"
REF2        = "9.3.cand02"
SIM_REF2    = "Yes"
RUN_STAT    = "Yes"
PHYSICS     = "FTFP_BERT"
EVENTS      = "10"
BFIELD      = "0"
#***endLOOKHERE***

# ---------------------------------------------

#***LOOKHERE***
CALORIMETER = "PbSci"
PARTICLE    = "pi+"
ENERGY      = "9GeV"
#***endLOOKHERE***

resultCode = os.system( "./simuDriver.sh " +
                        REF1 + " " + SIM_REF1 + " " +
                        REF2 + " " + SIM_REF2 + " " +
                        RUN_STAT + " " +
                        PHYSICS + " " +
                        CALORIMETER + " " +
                        PARTICLE + " " + ENERGY + " " +
                        EVENTS + " " + BFIELD )

if ( resultCode != 0 ) :
    print ' ***ERROR*** from: os.system( ./simuDriver.sh ... ) ! code=', resultCode
    
# ---------------------------------------------

print '========== END mainScript.py =========='
