# $Id: GNUmakefile,v 1.1 2000-01-07 14:48:45 grichine Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := TestEm8
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

LDFLAGS  += -L/cern/pro/lib
LOADLIBS += -lpacklib $(FCLIBS)

include $(G4INSTALL)/config/binmake.gmk
