# $Id: GNUmakefile,v 1.2 2000-01-19 20:55:07 vnivanch Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# G4NOHIST variable is checked.  Vladimir Ivanchenko, 2000/01/18
# --------------------------------------------------------------

name := TestEm6
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

ifdef G4NOHIST
  CPPFLAGS += -DG4NOHIST
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

ifndef G4NOHIST
  LDFLAGS  += -L/cern/pro/lib
  LOADLIBS += -lpacklib $(FCLIBS)
endif

include $(G4INSTALL)/config/binmake.gmk
