# $Id: GNUmakefile,v 1.7 2010-01-12 09:42:21 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for tests module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := test47
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../..
endif

.PHONY: all 

all: lib bin

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
LDFLAGS  += $(shell $(ROOTSYS)/bin/root-config --glibs)

include $(G4INSTALL)/config/binmake.gmk

