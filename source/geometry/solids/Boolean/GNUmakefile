# $Id: GNUmakefile,v 1.1 1999-01-07 16:07:48 gunter Exp $
# ----------------------------------------------------------------
# GNUmakefile for geometry/CSG library.  Gabriele Cosmo, 16/11/96.
# ----------------------------------------------------------------

name := G4geomBoolean

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include

include $(G4INSTALL)/config/common.gmk
