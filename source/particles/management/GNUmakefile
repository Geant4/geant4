# $Id: GNUmakefile,v 1.4 2000-10-24 13:36:51 gcosmo Exp $
# ----------------------------------------------------------------
# GNUmakefile for particle management library.  G.Cosmo, 4/4/97.
# ----------------------------------------------------------------

name := G4partman

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/materials/include


include $(G4INSTALL)/config/common.gmk
