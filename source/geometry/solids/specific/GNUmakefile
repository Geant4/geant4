# $Id: GNUmakefile,v 1.2 2000-04-27 10:37:26 gcosmo Exp $
# ----------------------------------------------------------------------------
# GNUmakefile for geometry/solids/specific library.  Gabriele Cosmo, 05/04/00.
# ----------------------------------------------------------------------------

name := G4specsolids

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include

include $(G4INSTALL)/config/common.gmk
