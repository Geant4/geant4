# $Id: GNUmakefile 75844 2013-11-06 17:35:56Z gcosmo $
# ------------------------------------------------------------
# GNUmakefile for phonons library. M.Kelsey 24-Oct-2013
# ------------------------------------------------------------

name := G4phonon

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/materials/include

include $(G4INSTALL)/config/common.gmk
