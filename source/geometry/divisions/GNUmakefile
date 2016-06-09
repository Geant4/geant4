# $Id: GNUmakefile,v 1.3 2003/11/19 11:51:44 gcosmo Exp $
# ----------------------------------------------------------------------
# GNUmakefile for geometry/divisions library.  Gabriele Cosmo, 16/06/03.
# ----------------------------------------------------------------------

name := G4geomdivision

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

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
