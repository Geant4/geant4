# --------------------------------------------------------------
# GNUmakefile for chips library. HPW 20-Nov-99
# --------------------------------------------------------------

name := G4hadronic_chips

SUBDIRS  =  body
SUBDIRS +=  nucleus
SUBDIRS +=  utils

SUBLIBS  = G4hadronic_body_ci
ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

include $(G4INSTALL)/config/globlib.gmk
