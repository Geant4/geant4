# $Id: GNUmakefile,v 1.3 2000-05-17 10:56:21 maire Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := TestEm4
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

ifdef G4NOHIST
  CPPFLAGS += -DG4NOHIST
else
  LDFLAGS  += -L$(CERN)/pro/lib
  LOADLIBS += -lpacklib $(FCLIBS)
endif

include $(G4INSTALL)/config/binmake.gmk
