# $Id: GNUmakefile,v 1.2 1999-07-13 14:19:42 graignac Exp $
# -----------------------------------------------------------------------
# GNUmakefile for geometry/management library.  Gabriele Cosmo, 16/11/96.
# -----------------------------------------------------------------------

name := G4geometrymng

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include 

#            -I$(G4BASE)/geometry/volumes/include

include $(G4INSTALL)/config/common.gmk
