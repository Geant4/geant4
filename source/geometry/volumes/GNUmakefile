# $Id: GNUmakefile,v 1.7 2004/06/11 14:17:19 gcosmo Exp $
# --------------------------------------------------------------------
# GNUmakefile for geometry/volumes library.  Gabriele Cosmo, 16/11/96.
# --------------------------------------------------------------------

name := G4volumes

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
            -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/geometry/verification/include \
            -I$(G4BASE)/geometry/management/include

ifdef G4DEBUG_NAVIGATION
  CPPFLAGS += -DG4DEBUG_NAVIGATION
endif
ifdef G4DEBUG_FIELD
  CPPFLAGS += -DG4DEBUG_FIELD
endif

include $(G4INSTALL)/config/common.gmk

