# $Id: GNUmakefile,v 1.1 2003/10/01 14:51:38 gcosmo Exp $
# --------------------------------------------------------------------
# GNUmakefile for geometry/volumes library.  Gabriele Cosmo, 16/11/96.
# --------------------------------------------------------------------

name := G4navigation

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/intercoms/include \
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

include $(G4INSTALL)/config/common.gmk

