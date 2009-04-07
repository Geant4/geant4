# $Id: GNUmakefile,v 1.4 2009-04-07 16:46:20 yarba Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := test47
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../
endif

.PHONY: all 

all: lib bin

#include $(G4INSTALL)/config/architecture.gmk
include misc/architecture.gmk

CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
LDFLAGS  += $(shell $(ROOTSYS)/bin/root-config --glibs)

#include $(G4INSTALL)/config/binmake.gmk
include misc/binmake.gmk

