# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.2 2002-06-13 09:25:36 francy Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------


G4TARGET := thyroid
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/binmake.gmk

# From lowenergy/test/GNUmakefile
CPPFLAGS += `aida-config --useHBook --include`
LDFLAGS  += `aida-config --useHBook --lib`
LOADLIBS  += `aida-config --useHBook --lib`














