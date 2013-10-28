# $Id: GNUmakefile 67085 2013-01-30 17:25:48Z yarba $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := CommonSW
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all 

all: lib 

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
LDFLAGS  += $(shell $(ROOTSYS)/bin/root-config --glibs)

#include $(G4INSTALL)/config/binmake.gmk
include ../../test48/binmake.gmk


