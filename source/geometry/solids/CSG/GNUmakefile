# $Id: GNUmakefile 101121 2016-11-07 09:18:01Z gcosmo $
# ----------------------------------------------------------------
# GNUmakefile for geometry/CSG library.  Gabriele Cosmo, 16/11/96.
# ----------------------------------------------------------------

name := G4csg

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4GEOM_ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPNumerics/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include

ifdef G4CSGDEBUG
  CPPFLAGS += -DG4CSGDEBUG
endif
ifdef G4BBOX_EXTENT
  CPPFLAGS += -DG4BBOX_EXTENT
endif

include $(G4INSTALL)/config/common.gmk
