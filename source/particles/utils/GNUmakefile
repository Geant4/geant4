# $Id: GNUmakefile,v 1.1 2003/09/21 19:38:49 kurasige Exp $
# ----------------------------------------------------------------
# GNUmakefile for particle management library.  G.Cosmo, 4/4/97.
# ----------------------------------------------------------------

name := G4partutils

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/particles/management/include


include $(G4INSTALL)/config/common.gmk
