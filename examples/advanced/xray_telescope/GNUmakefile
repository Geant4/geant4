# --------------------------------------------------------------
# $Id: GNUmakefile,v 1.6 2001-12-10 17:39:27 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := XrayTel
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

G4VIS_USE = 1

include $(G4INSTALL)/config/architecture.gmk

# Copy the XrayTelAnalysis files according to which version
# is selected.
# As the dependency is checked before "all" is reached, we need
# to do some tricks here ...

ifdef G4ANALYSIS_AIDA_CONFIG_CFLAGS
# AIDA (version 2.2 for now) selected
  FOO := $(shell cp -p include/XrayTelAnalysis.hh-AIDA include/XrayTelAnalysis.hh )
  FOO := $(shell cp -p src/XrayTelAnalysis.cc-AIDA src/XrayTelAnalysis.cc )
else
# Anaphe selected 
  include anaphe.gmk
  FOO := $(shell cp -p include/XrayTelAnalysis.hh-Anaphe include/XrayTelAnalysis.hh )
  FOO := $(shell cp -p src/XrayTelAnalysis.cc-Anaphe src/XrayTelAnalysis.cc )
endif

# End Analysis.

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk

