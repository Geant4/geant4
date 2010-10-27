# $Id: GNUmakefile,v 1.4 2010-10-27 07:47:05 gcosmo Exp $
# ----------------------------------------------------------------
# GNUmakefile for particle management library.  G.Cosmo, 4/4/97.
# ----------------------------------------------------------------

name := G4partutils

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4PARTICLES_ALLOC_EXPORT
CPPFLAGS += \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/particles/management/include


include $(G4INSTALL)/config/common.gmk
