# $Id: GNUmakefile,v 1.1 1999-01-07 16:07:23 gunter Exp $
# ------------------------------------------------------------
# GNUmakefile for BREPS library.  Gabriele Cosmo, 15/11/96.
# ------------------------------------------------------------

name := G4brep

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS  += \
             -I$(G4BASE)/geometry/solids/STEP/include \
             -I$(G4BASE)/geometry/management/include \
             -I$(G4BASE)/geometry/volumes/include \
             -I$(G4BASE)/global/management/include \
             -I$(G4BASE)/global/HEPGeometry/include \
             -I$(G4BASE)/graphics_reps/include

include $(G4INSTALL)/config/common.gmk

