# $Id: GNUmakefile,v 1.3 2001-02-07 17:30:58 gcosmo Exp $
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

$(G4INSTALL)/tests/tools/lib/$(G4SYSTEM)/libG4TestTool.a:
	cd $(G4INSTALL)/tests/tools/lib; $(MAKE)

