# $Id: GNUmakefile,v 1.5 2003-11-12 18:34:03 vnivanch Exp $
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
LDFLAGS    += -L/home/vnivanch/Anaphe/specific/redhat73/gcc-2.95.2/CLHEP/1.8.0.2/lib
LDFLAGS    += `aida-config --lib`

include $(G4INSTALL)/config/binmake.gmk
