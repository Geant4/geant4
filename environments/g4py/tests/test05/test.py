#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#
# ==================================================================
import test05

print "AClass(2,3.)"
a= test05.AClass(2,3.)
print "a.i=", a.i

print ""
print "func2(10)=", test05.func2(10)
print "func2(10,100)=", test05.func2(10, 100)

print ""
print "func3(1)=", test05.func3(1)
print "func3(1.)=", test05.func3(1.)

print ""
print "AMethod()=", a.AMethod()
print "AMethod(1)=", a.AMethod(1)
print "AMethod(1,2.)=", a.AMethod(1,2.)

print ""
print "BMethod()=", a.BMethod()
print "BMethod(1.)=", a.BMethod(1.)

print ""
print "CMethod(1)=", a.CMethod(1)
print "CMethod(1,10.)=", a.CMethod(1,10.)
print "CMethod(1,10.,100.)=", a.CMethod(1,10.,100.)

print ""
try:
  a.AMethod(1,2,3)  # error
except:
  print "*** ERRORS!! ***"
