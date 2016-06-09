# $Id: GNUmakefile,v 1.3 2002-01-10 15:52:07 gcosmo Exp $
# --------------------------------------------------------------------
# GNUmakefile for geometry/magneticfield library.  G.Cosmo, 22/01/97.

name := G4magneticfield

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include

ifdef G4DEBUG_FIELD
  CPPFLAGS += -DG4DEBUG_FIELD
endif

include $(G4INSTALL)/config/common.gmk
