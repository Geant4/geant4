#!/usr/bin/python

#----------------------------------------------------------------
# This Python
#----------------------------------------------------------------

import os
import sys

generalCase = "LHEP-FeSci-p-20GeV-5k"
caseA = "ref01"
caseB = "ref03"

os.system( "python2.2 drivePlot.py " + caseA + " " + caseB + " " + generalCase )

