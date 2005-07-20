# $Id: GNUmakefile,v 1.5 2005-07-20 20:37:35 allison Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := test19
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

CPPFLAGS += -I$(G4BASE)/visualization/XXX/include
EXTRALIBS += -lG4visXXX

# Temporary conditional use.  This will go in config/G4VIS_USE.gmk...
ifdef G4VIS_USE_RAYTRACERX
  CPPFLAGS += -DG4VIS_USE_RAYTRACERX
endif

include $(G4INSTALL)/config/binmake.gmk
