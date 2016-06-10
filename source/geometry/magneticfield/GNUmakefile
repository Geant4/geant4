# $Id: GNUmakefile 69786 2013-05-15 09:38:51Z gcosmo $
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
