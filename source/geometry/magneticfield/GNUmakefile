# $Id: GNUmakefile,v 1.1 1999-01-07 16:07:00 gunter Exp $
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
            -I$(G4BASE)/geometry/solids/CSG/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include

include $(G4INSTALL)/config/common.gmk
