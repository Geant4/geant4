# $Id: GNUmakefile,v 1.5 2002-12-06 15:19:40 stesting Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 27/06/98.
# --------------------------------------------------------------

name := test33
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

CPPFLAGS  += -I$(G4INSTALL)/tests/tools/include
EXTRA_LINK_DEPENDENCIES := \
  $(G4INSTALL)/tests/tools/lib/$(G4SYSTEM)/libG4TestTool.a 

$(G4INSTALL)/tests/tools/lib/$(G4SYSTEM)/libG4TestTool.a:
	cd $(G4INSTALL)/tests/tools/lib; $(MAKE)

.PHONY: all

all: lib bin

include $(G4INSTALL)/config/binmake.gmk

