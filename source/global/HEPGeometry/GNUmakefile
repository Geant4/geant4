# $Id: GNUmakefile,v 1.2 2004/06/11 14:18:56 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for global HEPGeometry.  Gabriele Cosmo, 26/9/96.
# --------------------------------------------------------------

name := G4hepgeometry

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4IOS_EXPORT
CPPFLAGS += -I$(G4BASE)/global/management/include

include $(G4INSTALL)/config/common.gmk
