#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#
# ==================================================================
import test04, sys

atrack= test04.Track()
print "<<< ref#(atrack)=", sys.getrefcount(atrack)


astep= atrack.GetStep1()
print "<<< ref#(atrack, astep)=",\
      sys.getrefcount(atrack),\
      sys.getrefcount(astep)
del astep
print "<<< delete astep"
print "<<< ref#(atrack)=", sys.getrefcount(atrack)

print ""
bstep= atrack.GetStep1()
bstep= test04.Step()
print "<<< previous step is still remain?"

print ""
print "--- stderr ---"

