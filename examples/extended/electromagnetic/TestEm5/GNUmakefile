# $Id: GNUmakefile,v 1.5 2001-02-21 14:44:32 maire Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := TestEm5
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

G4NOHIST := true
ifdef G4NOHIST
  CPPFLAGS += -DG4NOHIST
else
  LDFLAGS  += -L$(CERN)/pro/lib
  LOADLIBS += -lpacklib $(FCLIBS)
endif

include $(G4INSTALL)/config/binmake.gmk
