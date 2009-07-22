# $Id: GNUmakefile,v 1.5 2009-07-22 22:36:59 banerjee Exp $
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

include $(G4INSTALL)/config/architecture.gmk
#include misc/architecture.gmk

CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
LDFLAGS  += $(shell $(ROOTSYS)/bin/root-config --glibs)

include $(G4INSTALL)/config/binmake.gmk
#include misc/binmake.gmk

