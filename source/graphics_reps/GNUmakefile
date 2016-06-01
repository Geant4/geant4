# $Id: GNUmakefile,v 2.2 1998/11/06 14:01:53 allison Exp $
# --------------------------------------------------------------
# GNUmakefile for graphics_reps library.  John Allison, 19/7/96.
# --------------------------------------------------------------

name := G4graphics_reps

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPGeometry/include

include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib

