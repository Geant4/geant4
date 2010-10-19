# $Id: GNUmakefile,v 1.7 2010-10-19 10:18:01 gcosmo Exp $
# ----------------------------------------------------------------------------
# GNUmakefile for geometry/solids/specific library.  Gabriele Cosmo, 05/04/00.
# ----------------------------------------------------------------------------

name := G4specsolids

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4ALLOC_EXPORT
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
ifdef G4TWISTDEBUG
  CPPFLAGS += -DG4TWISTDEBUG
endif

include $(G4INSTALL)/config/common.gmk
