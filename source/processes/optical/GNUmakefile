# $Id: GNUmakefile,v 1.3 2004-10-28 14:26:38 gcosmo Exp $
# ------------------------------------------------------------
# GNUmakefile for optical photons library. G.Folger 10-Dec-97.
# ------------------------------------------------------------

name := G4optical

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/navigation/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/materials/include

ifdef G4OPTICAL_DEBUG
  CPPFLAGS += -DG4OPTICAL_DEBUG
endif

include $(G4INSTALL)/config/common.gmk
