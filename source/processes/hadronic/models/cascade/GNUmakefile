# --------------------------------------------------------------
# GNUmakefile for chips library. HPW 20-Nov-99
# --------------------------------------------------------------

name := G4hadronic_hetcpp

SUBDIRS  =  evaporation
SUBDIRS  +=  utils

SUBLIBS  = G4hadronic_hetcpp_evaporation
SUBLIBS  += G4hadronic_hetcpp_utils

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif


include $(G4INSTALL)/config/globlib.gmk
