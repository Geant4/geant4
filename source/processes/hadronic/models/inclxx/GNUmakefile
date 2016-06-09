# --------------------------------------------------------------
# GNUmakefile for chips library. HPW 20-Nov-99
# --------------------------------------------------------------

name := G4had_inclxx

SUBDIRS  = utils
SUBDIRS  += incl_physics 
SUBDIRS  += interface

SUBLIBS  = G4hadronic_inclxx_utils
SUBLIBS  += G4hadronic_inclxx_physics
SUBLIBS  += G4hadronic_inclxx_interface

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif


include $(G4INSTALL)/config/globlib.gmk
