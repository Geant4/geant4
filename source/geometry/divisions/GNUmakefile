# $Id: GNUmakefile 79490 2014-03-05 15:23:33Z gcosmo $
# ----------------------------------------------------------------------
# GNUmakefile for geometry/divisions library.  Gabriele Cosmo, 16/06/03.
# ----------------------------------------------------------------------

name := G4geomdivision

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4GEOM_ALLOC_EXPORT
CPPFLAGS += \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/geometry/solids/CSG/include \
            -I$(G4BASE)/geometry/solids/specific/include \
            -I$(G4BASE)/geometry/management/include

ifdef G4DIVDEBUG
  CPPFLAGS += -DG4DIVDEBUG
endif

include  $(G4INSTALL)/config/common.gmk
