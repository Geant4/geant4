# $Id: GNUmakefile,v 1.10 2001-12-07 07:45:42 pfeiffer Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := GammaRayTel
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: tmpana lib bin

include $(G4INSTALL)/config/architecture.gmk


# Check whether to use AIDA (default) or not.
# Copy the GammaRayTelAnalysis files according to which version
# is selected. For this the target tmpana has been created.
# This part of the makefile will be obsolete soon (see also below).

ifdef G4ANALYSIS_AIDA_CONFIG_CFLAGS
# AIDA (version 2.2 for now) selected
tmpana:
	cp include/GammaRayTelAnalysis.hh-AIDA include/GammaRayTelAnalysis.hh
	cp src/GammaRayTelAnalysis.cc-AIDA src/GammaRayTelAnalysis.cc
else
# Anaphe selected (will be obsolete as soon as Anaphe is AIDA 2.2 compliant)
include anaphe.gmk
tmpana:
	cp include/GammaRayTelAnalysis.hh-Anaphe include/GammaRayTelAnalysis.hh
	cp src/GammaRayTelAnalysis.cc-Anaphe src/GammaRayTelAnalysis.cc
endif

# End check for AIDA.


include $(G4INSTALL)/config/binmake.gmk





