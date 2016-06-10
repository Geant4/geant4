# $Id: GNUmakefile 66892 2013-01-17 10:57:59Z gunter $
# --------------------------------------------------------------
# GNUmakefile for global management.  Gabriele Cosmo, 26/9/96.
# --------------------------------------------------------------

name := G4globman

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4GLOB_ALLOC_EXPORT
ifdef G4RF_DEBUG
  CPPFLAGS += -DG4RF_DEBUG
endif

include $(G4INSTALL)/config/common.gmk
