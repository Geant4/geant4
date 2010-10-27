# $Id: GNUmakefile,v 1.6 2010-10-27 07:34:32 gcosmo Exp $
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
            -I$(G4BASE)/geometry/solids/CSG/include \
            -I$(G4BASE)/geometry/solids/specific/include \
            -I$(G4BASE)/geometry/management/include

ifdef G4DIVDEBUG
  CPPFLAGS += -DG4DIVDEBUG
endif

include  $(G4INSTALL)/config/common.gmk
