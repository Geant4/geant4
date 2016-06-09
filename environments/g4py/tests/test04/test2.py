#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#
# ==================================================================
import test04, sys

btrack= test04.Track()
print "<<< ref#(btrack)=", sys.getrefcount(btrack)
bstep= btrack.GetStep2()
print "<<< ref#(btrack, bstep)=",\
      sys.getrefcount(btrack),\
      sys.getrefcount(bstep)
del bstep
print "<<< delete bstep"
print "<<< ref#(btrack)=", sys.getrefcount(btrack)


print ""
cstep= btrack.GetStep2()
cstep= test04.Step()
print "<<< previous step is still remain?"


print ""
print "--- stderr ---"


