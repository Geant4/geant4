# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.3 2002-06-13 17:15:22 francy Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name :=thyroid
G4TARGET := $(name)
G4EXLIB := true
G4ANALYSIS_USE=1
ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

#include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/binmake.gmk

# From lowenergy/test/GNUmakefile
#CPPFLAGS += `aida-config --useHBook --include`
#LDFLAGS  += `aida-config --useHBook --lib`
#LOADLIBS  += `aida-config --useHBook --lib`
#LOADLIBS  += `aida-config --useHBook --lib`













