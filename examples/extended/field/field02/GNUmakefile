# $Id: GNUmakefile,v 1.2 2001-03-27 16:25:03 grichine Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := field02
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
  LDFLAGS  += -L/cern/pro/lib
  LOADLIBS += -lpacklib $(FCLIBS)
endif


include $(G4INSTALL)/config/binmake.gmk
