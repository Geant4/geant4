# $Id: GNUmakefile 66356 2012-12-18 09:02:32Z gcosmo $
# --------------------------------------------------------------------
# GNUmakefile for geometry/usolids library.  Gabriele Cosmo, 16/11/96.
# --------------------------------------------------------------------

name := G4geomUSolids

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4GEOM_ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/management/include \

ifdef G4USOLDEBUG
  CPPFLAGS += -DG4USOLDEBUG
endif

include $(G4INSTALL)/config/common.gmk
