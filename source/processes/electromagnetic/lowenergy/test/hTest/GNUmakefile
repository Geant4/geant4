# $Id: GNUmakefile,v 1.6 2002-12-18 14:29:27 vnivanch Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# G4NOHIST variable is checked.  Vladimir Ivanchenko, 2000/01/18
# --------------------------------------------------------------

name := hTest
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../../../..
endif

.PHONY: all
all: lib bin

ifndef G4ANALYSIS_USE
  G4ANALYSIS_USE = 1
endif

CPPFLAGS += -DG4ANALYSIS_USE

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += `aida-config --include`
LDFLAGS  += `aida-config --lib`

include $(G4INSTALL)/config/binmake.gmk




