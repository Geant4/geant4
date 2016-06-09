# $Id:$
# ------------------------------------------------------------
# GNUmakefile for internal CLHEP library.  Gabriele Cosmo, 4/5/11.
# ------------------------------------------------------------

name := G4clhep

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

GLOBLIBS  = 

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DCLHEP_EXPORT

include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib

