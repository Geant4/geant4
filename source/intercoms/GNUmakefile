# $Id: GNUmakefile 106172 2017-09-15 13:03:57Z gcosmo $
# -------------------------------------------------------------
# GNUmakefile for intercoms library.  Gabriele Cosmo, 11/11/96.
# -------------------------------------------------------------

name := G4intercoms

ifndef G4INSTALL
  G4INSTALL = ../..
endif

GLOBLIBS = libG4global.lib

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4ICOMS_ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include

include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib

