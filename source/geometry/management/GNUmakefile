# $Id: GNUmakefile,v 1.5 2004/11/12 09:05:05 gcosmo Exp $
# -----------------------------------------------------------------------
# GNUmakefile for geometry/management library.  Gabriele Cosmo, 16/11/96.
# -----------------------------------------------------------------------

name := G4geometrymng

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4GEOM_ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/materials/include

ifdef G4GEOMETRY_VOXELDEBUG
  CPPFLAGS += -DG4GEOMETRY_VOXELDEBUG
endif

include $(G4INSTALL)/config/common.gmk
