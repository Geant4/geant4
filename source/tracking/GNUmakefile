# $Id: GNUmakefile,v 1.3 2003/10/03 14:50:21 gcosmo Exp $
# ----------------------------------------------------------
# GNUmakefile for tracking library.  Katsuya Amako, 5/9/95.
# ----------------------------------------------------------

name := G4tracking

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/navigation/include \
            -I$(G4BASE)/geometry/solids/CSG/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/cuts/include \
            -I$(G4BASE)/processes/electromagnetic/utils/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/digits_hits/detector/include \
            -I$(G4BASE)/digits_hits/hits/include \
            -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/intercoms/include

include $(G4INSTALL)/config/common.gmk

.PHONY: global

global: lib
