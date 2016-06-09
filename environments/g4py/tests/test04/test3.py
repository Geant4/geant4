#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
# This script is just for test use. It will crash.
# ==================================================================
import test04, sys

print "*** This script is just for test use. It will crash!"

ctrack= test04.Track()
print "<<< ref#(ctrack)=", sys.getrefcount(ctrack)
cstep= ctrack.GetStep3()
print "<<< ref#(ctrack, cstep)=",\
      sys.getrefcount(ctrack),\
      sys.getrefcount(cstep)
del cstep
print "<<< delete cstep"
print "<<< ref#(ctrack)=", sys.getrefcount(ctrack)


print ""
cstep= ctrack.GetStep3()
cstep= test04.Step()
print "<<< previous step is still remain?"


print ""
print "--- stderr ---"


