# $Id: GNUmakefile,v 1.1 1999-01-08 16:34:27 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------

name := test01
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

CPPFLAGS+=-I$(G4INSTALL)/tests/tools/include
EXTRALIBS+=-L$(G4INSTALL)/tests/tools/lib/$(G4SYSTEM) -lG4TestTool

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
