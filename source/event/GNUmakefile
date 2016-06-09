# $Id: GNUmakefile,v 1.5 2003/10/03 14:12:14 gcosmo Exp $
# ------------------------------------------------------------
# GNUmakefile for events library.  Makoto Asai, 5/9/95.
# ------------------------------------------------------------

name := G4event

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/global/HEPNumerics/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/tracking/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/electromagnetic/utils/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/navigation/include \
	    -I$(G4BASE)/geometry/biasing/include \
	    -I$(G4BASE)/digits_hits/detector/include \
	    -I$(G4BASE)/digits_hits/hits/include \
	    -I$(G4BASE)/digits_hits/digits/include \
	    -I$(G4BASE)/intercoms/include

ifdef G4_STORE_TRAJECTORY
  CPPFLAGS  += -DG4_STORE_TRAJECTORY
endif

include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib

