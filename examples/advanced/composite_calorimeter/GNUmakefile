# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.4 2002-10-30 20:00:16 ribon Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := CompositeCalorimeter

G4TARGET := $(name)
G4EXLIB := true

G4ANALYSIS_USE = 1

ifndef G4INSTALL
  G4INSTALL = ../../../
endif

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

CPPFLAGS += `aida-config --include`
LDFLAGS  += `aida-config --lib`
LOADLIBS  += `aida-config --lib`


