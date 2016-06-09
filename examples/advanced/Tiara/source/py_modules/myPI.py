# $Id: myPI.py,v 1.1 2004/06/09 15:04:36 daquinog Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-09-02 $
# -------------------------------------------------------------------
#
import sys

origArgv = sys.argv
sys.argv = ["","--noGraphics"]
from PyAida import *
sys.argv = origArgv
