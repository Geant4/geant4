# $Id: GNUmakefile 101118 2016-11-07 09:10:59Z gcosmo $
# ----------------------------------------------------------------------------
# GNUmakefile for geometry/solids/specific library.  Gabriele Cosmo, 05/04/00.
# ----------------------------------------------------------------------------

name := G4specsolids

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4GEOM_ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPNumerics/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/solids/CSG/include \
            -I$(G4BASE)/geometry/management/include

ifdef G4SPECSDEBUG
  CPPFLAGS += -DG4SPECSDEBUG
endif
ifdef G4TWISTDEBUG
  CPPFLAGS += -DG4TWISTDEBUG
endif
ifdef G4BBOX_EXTENT
  CPPFLAGS += -DG4BBOX_EXTENT
endif

include $(G4INSTALL)/config/common.gmk
