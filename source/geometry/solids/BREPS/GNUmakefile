# $Id: GNUmakefile,v 1.6 2010-10-19 10:18:01 gcosmo Exp $
# ------------------------------------------------------------
# GNUmakefile for BREPS library.  Gabriele Cosmo, 15/11/96.
# ------------------------------------------------------------

name := G4brep

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4ALLOC_EXPORT
CPPFLAGS += \
             -I$(G4BASE)/geometry/management/include \
             -I$(G4BASE)/geometry/solids/Boolean/include \
             -I$(G4BASE)/geometry/solids/CSG/include \
             -I$(G4BASE)/global/management/include \
             -I$(G4BASE)/global/HEPGeometry/include \
             -I$(G4BASE)/graphics_reps/include \
             -I$(G4BASE)/intercoms/include

include $(G4INSTALL)/config/common.gmk
