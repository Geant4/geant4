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
#               -  the reference to be run: eg:   "9.1.p02"
#               -  the Physics Lists:
#                  "LHEP", "QGSP", "QGSC", "QGSP_BERT", "QGSP_BIC", etc.
#               -  the Number of Events:
#                  e.g. : "5000", "5k", etc.
#               -  the Bfield:
#                  e.g. : "0", "4tesla", etc.
#
REF         = "9.1.p02"
PHYSICS     = "QGSP"
EVENTS      = "10"
BFIELD      = "4tesla"
#***endLOOKHERE***

# ---------------------------------------------

resultCode = os.system( "./simuDriver.sh " +
                        REF + " " +
                        PHYSICS + " " +
                        EVENTS + " " + BFIELD )

if ( resultCode != 0 ) :
    print ' ***ERROR*** from: os.system( ./simuDriver.sh ... ) ! code=', resultCode
    
# ---------------------------------------------

print '========== END mainScript.py =========='
