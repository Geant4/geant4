# $Id: GNUmakefile 99159 2016-09-07 08:11:50Z gcosmo $
# ----------------------------------------------------------------
# GNUmakefile for particle management library.  G.Cosmo, 4/4/97.
# ----------------------------------------------------------------

name := G4partman

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
            -I$(G4BASE)/materials/include


include $(G4INSTALL)/config/common.gmk

ifdef G4NucleiProperties_USE_OLD_AME_TABLE
   CPPFLAGS += -DG4NucleiProperties_USE_OLD_AME_TABLE
endif
