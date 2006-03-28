# $Id: GNUmakefile,v 1.7 2006-03-28 15:36:27 allison Exp $
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

ifdef G4VIS_USE_XXXSG
  CPPFLAGS += -DG4VIS_USE_XXXSG
endif

include $(G4INSTALL)/config/binmake.gmk
