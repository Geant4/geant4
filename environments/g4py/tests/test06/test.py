#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#
# ==================================================================
import test06

class MyZClass1(test06.ZBase):
  "My Class derived from ZBase(C++)"
  def AMethod(self):
    print "MyZClass::AMethod is called:"


class MyZClass2(test06.ZBase):
  "My Class derived from ZBase(C++)"
  def VMethod(self, message):
    print "MyZClass::VMethod is called:", message


class MyXClass(test06.XBase):
  "My Class derived from XBase(C++)"
  def PVMethod(self):
    return "PVMethod is called from MyClass"


print "***  ZBase..."
z= test06.ZBase()
z.AMethod()
z.VMethod("hello world")

print ""
print "*** Inherit from ZBase class..."
myz1= MyZClass1()
myz1.AMethod()
myz1.VMethod("XXXX")

print "*** w/ override"
myz2= MyZClass2()
myz2.VMethod("YYYY")


print ""
print "*** Inherit from XBase, pure virtual class..."
myx= MyXClass()
print "myx.ival=", myx.ival
print myx.PVMethod()

print ""
x= test06.XBase()
print "*** Runtime ERROR will occur!" \
      " because pure virtual function is called."
x.PVMethod()


