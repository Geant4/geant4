#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#
# ==================================================================
import test09


a= test09.AClass(10)
b= test09.AClass(20)

print "*** (a, b)=", a, b

c= a+b
print "*** a+b=", c

a+=b
print "*** a & b are merged, a=", a

print "*** a==b ?", (a==b)

