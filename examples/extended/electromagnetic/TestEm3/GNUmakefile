# $Id: GNUmakefile,v 1.1 1999-10-11 16:55:33 maire Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := TestEm3
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

LDFLAGS  += -L/cern/pro/lib
LOADLIBS += -lpacklib $(FCLIBS)

include $(G4INSTALL)/config/binmake.gmk
