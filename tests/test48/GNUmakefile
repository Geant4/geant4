# $Id$
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := test48
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../
endif

.PHONY: all 

all: lib bin

include $(G4INSTALL)/config/architecture.gmk

ifeq (,$(ROOTSYS))
  $(error test48 cannot be built without ROOT.  Ensure that ROOTSYS is defined)
else
  CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
  LDFLAGS  += $(shell $(ROOTSYS)/bin/root-config --glibs)
endif

#include $(G4INSTALL)/config/binmake.gmk
include binmake.gmk

