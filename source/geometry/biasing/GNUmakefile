# $Id: GNUmakefile,v 1.1 2002-03-25 17:09:16 gcosmo Exp $
# --------------------------------------------------------------------
# GNUmakefile for geometry/biasing library.  Gabriele Cosmo, 25/03/02.
# --------------------------------------------------------------------

name := G4geombias

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include

include $(G4INSTALL)/config/common.gmk
