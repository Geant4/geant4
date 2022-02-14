# ------------------------------------------------------------
# GNUmakefile for internal PTL library.  Gabriele Cosmo, 18/6/20.
# ------------------------------------------------------------

name := G4ptl

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

GLOBLIBS  = 

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DPTL_ALLOC_EXPORT

include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib

