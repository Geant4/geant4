# $Id: GNUmakefile,v 1.1 1999-01-08 16:35:06 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := test07
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all
all: lib bin

ifdef MAKEHBOOK
 CPPFLAGS+=-I$(G4INSTALL)/tests/tools/include -DMAKEHISTO
 EXTRALIBS+=-L$(G4INSTALL)/tests/tools/lib/$(G4SYSTEM) -lG4TestTool
 EXTRALIBS+=-L/cern/pro/lib -lpacklib
ifeq ($(G4SYSTEM),G++-Linux)
 EXTRALIBS+=-lf2c -lgcc
endif
ifeq ($(G4SYSTEM),AIX-AFS)
 EXTRALIBS+=-lxlf90 -lxlf
endif
else
CPPFLAGS+=-I$(G4INSTALL)/tests/tools/include
EXTRALIBS+=-L$(G4INSTALL)/tests/tools/lib/$(G4SYSTEM) -lG4TestTool
endif

include $(G4INSTALL)/config/binmake.gmk

