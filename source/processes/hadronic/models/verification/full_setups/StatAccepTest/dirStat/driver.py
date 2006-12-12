#!/usr/bin/python

#---------------------------------------------------------------------
# This Python script has 3 input parameters:
#   1) a Geant4 reference tag ;       example:  4.6.2.ref03
#   2) another Geant4 reference tag ; example:  4.6.2.ref04
#   3) a string summarizing which Physics List,
#      which calorimeter setup,
#      which beam particle, which beam energy,
#      which  number of events;
#                                     example:  LHEP-FeSci-p-20GeV-5k
#
# This script then drives the Python script  drivePlot.py .
#---------------------------------------------------------------------

import os
import sys

print '   ========== BEGIN driver.py ========== '

caseA = sys.argv[1]
caseB = sys.argv[2]
generalCase = sys.argv[3]

resultCode = os.system( "python drivePlot.py " +
                        caseA + " " + caseB + " " + generalCase )

if ( resultCode != 0 ) :
    print ' ***ERROR*** from: os.system( python drivePlot.py ... ) ! code=', resultCode
    sys.exit( 41 )

print '   ========== END driver.py ========== '
