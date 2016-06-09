#!/usr/bin/python
# ==================================================================
# python script for Geant4Py test
#
#
# ==================================================================
import test07

acolor= test07.Color().red

print acolor==test07.Color().red
print acolor==test07.Color().blue
print acolor==test07.Color().yellow

