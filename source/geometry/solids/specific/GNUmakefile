# $Id: GNUmakefile,v 1.5 2004/12/17 16:34:55 link Exp $
# ----------------------------------------------------------------------------
# GNUmakefile for geometry/solids/specific library.  Gabriele Cosmo, 05/04/00.
# ----------------------------------------------------------------------------

name := G4specsolids

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4GEOM_ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPNumerics/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include

ifdef G4SPECSDEBUG
  CPPFLAGS += -DG4SPECSDEBUG
endif

include $(G4INSTALL)/config/common.gmk
