# $Id: GNUmakefile,v 1.1 2002-03-18 17:32:48 vnivanch Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := test30
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../
endif

.PHONY: all tools
 

all: lib bin

include $(G4INSTALL)/config/architecture.gmk

ifdef G4NOHIST
  CPPFLAGS += -DG4NOHIST
else
  LDFLAGS  += -L/cern/pro/lib
  LOADLIBS += -lpacklib $(FCLIBS)
endif


include $(G4INSTALL)/config/binmake.gmk


