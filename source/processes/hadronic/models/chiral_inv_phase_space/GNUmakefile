# --------------------------------------------------------------
# GNUmakefile for chips library. HPW 20-Nov-99
# --------------------------------------------------------------

name := G4hadronic_chips

SUBDIRS  =  body
SUBDIRS  =  interface

SUBLIBS  = G4hadronic_body_ci
SUBLIBS  = G4hadronic_interface_ci

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

include $(G4INSTALL)/config/globlib.gmk
