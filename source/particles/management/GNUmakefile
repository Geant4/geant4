# $Id: GNUmakefile,v 1.2 1999-10-29 08:03:48 kurasige Exp $
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

# ----------
# for using STL map for G4ParticlTable
#
ifdef G4USE_STL
 CPPFLAGS += -DG4USE_STL_MAP
endif 

include $(G4INSTALL)/config/common.gmk
