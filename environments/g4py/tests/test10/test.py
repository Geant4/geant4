#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#
# ==================================================================
import test10

class MyX(test10.XBase):
  "My X"

  def VMethodA(self, a):
    print "*** MyX::VMethod...A() is called."
    
  def VMethodB(self, b):
    print "*** MyX::VMethod...B() is called."


x= test10.XBase()
myx= MyX()
z= test10.ZClass()

print "!!! direct call via C++ pointer"
z.SetXBase(x)
z.Process()

print ""
print "!!! call via python inheritor"
z.SetXBase(myx)
z.Process()

