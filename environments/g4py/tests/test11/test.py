#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#
# ==================================================================
import test11

a= test11.AClass(1)

try:
  b = test11.AClass()
except:
  print "AClass default constructor is not allowed."
