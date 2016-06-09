# $Id: GNUmakefile,v 1.5.10.1 2009/03/03 13:58:24 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for global management.  Gabriele Cosmo, 26/9/96.
# --------------------------------------------------------------

name := G4globman

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4IOS_EXPORT
ifdef G4RF_DEBUG
  CPPFLAGS += -DG4RF_DEBUG
endif

include $(G4INSTALL)/config/common.gmk
