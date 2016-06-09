# --------------------------------------------------------------
# GNUmakefile for chips library. HPW 20-Nov-99
# --------------------------------------------------------------

name := G4had_hetcpp

SUBDIRS  =  evaporation
SUBDIRS  +=  utils
SUBDIRS  +=  cascade

SUBLIBS  = G4hadronic_hetcpp_evaporation
SUBLIBS  += G4hadronic_hetcpp_utils
SUBLIBS  += G4hadronic_bert_cascade

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif


include $(G4INSTALL)/config/globlib.gmk
