# $Id: GNUmakefile,v 1.5 2002-07-24 16:45:34 vnivanch Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# G4NOHIST variable is checked.  Vladimir Ivanchenko, 2000/01/18
# --------------------------------------------------------------

name := hTest
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

#
# NEW: this requires the env-var LHCXX_REL_DIR to be set up correctly,
# see the first lines in G4ComptonTest.cc
#
OGLLIBS   := -L$(OGLHOME)/lib -lGLU -lGL
LOADLIBS += $(OGLLIBS)
CPPFLAGS += `aida-config --useHBook --include`
LDFLAGS  += `aida-config --useHBook --lib`
#EXTRALIBS += -L$(OGLHOME)/lib -lGL

include $(G4INSTALL)/config/binmake.gmk




