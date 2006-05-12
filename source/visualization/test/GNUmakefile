# $Id: GNUmakefile,v 1.8 2006-05-12 12:39:16 allison Exp $
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

include $(G4INSTALL)/config/binmake.gmk

ifdef G4VIS_USE
CPPFLAGS += -I$(G4BASE)/visualization/XXX/include
EXTRALIBS += -lG4visXXX
endif

ifdef G4VIS_USE_XXXSG
  CPPFLAGS += -DG4VIS_USE_XXXSG
endif
