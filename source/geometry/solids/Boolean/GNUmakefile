# $Id: GNUmakefile,v 1.5 2005-03-04 09:26:54 grichine Exp $
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
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include

ifdef G4BOOLDEBUG
  CPPFLAGS += -DG4BOOLDEBUG
endif
ifdef G4NEW_SURF_NORMAL
  CPPFLAGS += -DG4NEW_SURF_NORMAL
endif

include $(G4INSTALL)/config/common.gmk
