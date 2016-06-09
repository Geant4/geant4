#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#
# ==================================================================
import test01

print "importing created module...\n"

a= test01.AClass()
b= test01.AClass(0)
c= test01.AClass(1, 1.)

print " >>a<<",
a.AMethod()

print ">>b<<",
b.AMethod()

print ">>c<<",
c.AMethod()

