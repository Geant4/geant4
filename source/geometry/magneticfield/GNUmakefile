# $Id: GNUmakefile,v 1.2 1999-02-17 16:57:35 japost Exp $
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
            -I$(G4BASE)/global/HEPGeometry/include \

include $(G4INSTALL)/config/common.gmk
