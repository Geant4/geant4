# $Id: GNUmakefile,v 1.3 2003-01-14 11:58:25 vnivanch Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := test31
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../
endif

.PHONY: all
all: lib bin

CPPFLAGS += -DG4ANALYSIS_USE

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS   += `aida-config --include`
LDFLAGS    += `aida-config --lib`

include $(G4INSTALL)/config/binmake.gmk