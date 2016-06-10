# $Id: GNUmakefile 66892 2013-01-17 10:57:59Z gunter $
# -------------------------------------------------------------
# GNUmakefile for intercoms library.  Gabriele Cosmo, 11/11/96.
# -------------------------------------------------------------

name := G4intercoms

ifndef G4INSTALL
  G4INSTALL = ../..
endif

GLOBLIBS = libG4global.lib

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include

include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib

