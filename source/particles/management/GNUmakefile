# $Id: GNUmakefile,v 1.1 1999-01-07 16:10:26 gunter Exp $
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
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/decay/include \
            -I$(G4BASE)/processes/electromagnetic/utils/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/materials/include

include $(G4INSTALL)/config/common.gmk
