# $Id: GNUmakefile,v 1.2 2001-02-01 21:24:23 radoone Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------

name := testAssemblyVolume
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

CPPFLAGS+=-I$(G4INSTALL)/tests/tools/include
EXTRALIBS+=$(LIB_PATH)$(G4INSTALL)/tests/tools/lib/$(G4SYSTEM) -lG4TestTool

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
