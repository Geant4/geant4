# $Id: GNUmakefile,v 1.6 2010-10-19 10:18:01 gcosmo Exp $
# --------------------------------------------------------------------
# GNUmakefile for geometry/volumes library.  Gabriele Cosmo, 16/11/96.
# --------------------------------------------------------------------

name := G4navigation

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include

ifdef G4DEBUG_NAVIGATION
  CPPFLAGS += -DG4DEBUG_NAVIGATION
endif
ifdef G4DEBUG_FIELD
  CPPFLAGS += -DG4DEBUG_FIELD
endif
ifdef G4DEBUG_PATHFINDER
  CPPFLAGS += -DG4DEBUG_PATHFINDER
endif
ifdef PATHFINDER_OPTIMISATION
  CPPFLAGS += -DPATHFINDER_OPTIMISATION
endif

include $(G4INSTALL)/config/common.gmk

