# $Id: myLiz.py,v 1.2 2003-06-16 17:06:44 dressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: not supported by cvs2svn $
# -------------------------------------------------------------------
#
import sys
origArgv = sys.argv
sys.argv = ["","--noGraphics"]
from liz import *
sys.argv = origArgv
