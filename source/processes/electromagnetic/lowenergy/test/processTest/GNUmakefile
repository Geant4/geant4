# $Id: GNUmakefile,v 1.1 2001-10-15 12:27:03 pia Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# G4NOHIST variable is checked.  Vladimir Ivanchenko, 2000/01/18
# --------------------------------------------------------------

name := processTest
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk

ifdef G4NOHIST
  CPPFLAGS += -DG4NOHIST
else
  LDFLAGS  += -L/cern/pro/lib
  LOADLIBS += -lpacklib $(FCLIBS)
endif

include $(G4INSTALL)/config/binmake.gmk




