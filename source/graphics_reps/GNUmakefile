# $Id: GNUmakefile,v 1.3 2005-03-26 22:37:14 allison Exp $
# --------------------------------------------------------------
# GNUmakefile for graphics_reps library.  John Allison, 19/7/96.
# --------------------------------------------------------------

name := G4graphics_reps

ifndef G4INSTALL
  G4INSTALL = ../..
endif

GLOBLIBS = libG4global.lib

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/global/HEPGeometry/include
CPPFLAGS += -I$(G4BASE)/intercoms/include

include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib

