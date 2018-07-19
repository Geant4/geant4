# $Id: GNUmakefile 103464 2017-04-11 07:24:03Z gcosmo $
# ----------------------------------------------------------------
# GNUmakefile for geometry/CSG library.  Gabriele Cosmo, 16/11/96.
# ----------------------------------------------------------------

name := G4geomBoolean

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4GEOM_ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/solids/CSG/include \
            -I$(G4BASE)/geometry/solids/specific/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include

ifdef G4BOOLDEBUG
  CPPFLAGS += -DG4BOOLDEBUG
endif

include $(G4INSTALL)/config/common.gmk
