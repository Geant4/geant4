# $Id: GNUmakefile,v 1.1 1999-01-07 16:09:08 gunter Exp $
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

