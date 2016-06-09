# $Id: GNUmakefile,v 1.2 2004/06/01 15:36:32 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for graphics_reps library.  John Allison, 19/7/96.
# --------------------------------------------------------------

name := G4graphics_reps

ifndef G4INSTALL
  G4INSTALL = ../..
endif

GLOBLIBS = libG4global.lib

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPGeometry/include

include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib

