# $Id: GNUmakefile,v 1.2 2002-10-29 15:43:02 dressel Exp $
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
EXTRALIBS += $(LIB_PATH)$(G4INSTALL)/tests/tools/lib/$(G4SYSTEM) -lG4TestTool -lXmu -lXm
EXTRA_LINK_DEPENDENCIES := \
  $(G4INSTALL)/tests/tools/lib/$(G4SYSTEM)/libG4TestTool.a 



.PHONY: all

all: lib bin

include $(G4INSTALL)/config/binmake.gmk

$(G4INSTALL)/tests/tools/lib/$(G4SYSTEM)/libG4TestTool.a:
	cd $(G4INSTALL)/tests/tools/lib; $(MAKE)
