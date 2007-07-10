#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#
# ==================================================================
import test12

a= test12.AClass()
a.ival=10

b= test12.AClass()
b.ival=20

c= test12.AClass()
c.ival=30

v= test12.AVector()
v[:]= [a, b, c]

print "*** size of vector= ", len(v)
x=v[0]

print ""
print "*** a vs v[0]"
a.Print()
x.Print()

print ""
print "*** dump vector..."
test12.PrintVector(v)





