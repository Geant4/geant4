# $Id: GNUmakefile,v 1.1 2003-06-16 15:11:40 gcosmo Exp $
# ----------------------------------------------------------------------
# GNUmakefile for geometry/divisions library.  Gabriele Cosmo, 16/06/03.
# ----------------------------------------------------------------------

name :=  G4geomdivision

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

include  $(G4INSTALL)/config/common.gmk
