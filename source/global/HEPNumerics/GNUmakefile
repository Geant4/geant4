# $Id: GNUmakefile,v 1.4 2010-10-27 07:40:06 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for global HEPNumerics.  Gabriele Cosmo, 26/9/96.
# --------------------------------------------------------------

name := G4hepnumerics

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4GLOB_ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/global/management/include

include $(G4INSTALL)/config/common.gmk
