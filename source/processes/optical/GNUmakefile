# $Id: GNUmakefile,v 1.1 1999-01-07 16:14:00 gunter Exp $
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
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/materials/include
           
include $(G4INSTALL)/config/common.gmk
