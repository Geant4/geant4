# $Id: GNUmakefile,v 1.4 2000-01-21 13:47:35 gcosmo Exp $
# ------------------------------------------------------------
# GNUmakefile for BREPS library.  Gabriele Cosmo, 15/11/96.
# ------------------------------------------------------------

name := G4brep

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS  += \
             -I$(G4BASE)/geometry/management/include \
             -I$(G4BASE)/geometry/solids/Boolean/include \
             -I$(G4BASE)/geometry/solids/CSG/include \
             -I$(G4BASE)/global/management/include \
             -I$(G4BASE)/global/HEPGeometry/include \
             -I$(G4BASE)/graphics_reps/include \
             -I$(G4BASE)/intercoms/include

include $(G4INSTALL)/config/common.gmk
