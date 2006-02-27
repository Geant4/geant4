#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#
# ==================================================================
import test02

print "importing created module...\n"

a= test02.AClass()
b= test02.BClass()

print ">>a<<",
a.AMethod()
print ""

print ">>b<<",
b.AMethod()
print ""

print "%%% a.VMethod(a)=", a.VMethod(a)
print "%%% a.VMethod(b)=", a.VMethod(b)

print "%%% b.VMethod(a)=", b.VMethod(a)
print "%%% b.VMethod(b)=", b.VMethod(b)

