# $Id: GNUmakefile,v 1.1 1999-01-07 16:14:11 gunter Exp $
# -----------------------------------------------------------
# GNUmakefile for digits+hits library.  Makoto Asai, 1/11/96.
# -----------------------------------------------------------

name := G4readout

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += \
	    -I$(G4BASE)/global/HEPRandom/include \
	    -I$(G4BASE)/global/management/include \
	    -I$(G4BASE)/global/HEPGeometry/include \
	    -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/electromagnetic/utils/include \
	    -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/tracking/include \
            -I$(G4BASE)/run/include \
            -I$(G4BASE)/digits+hits/detector/include \
            -I$(G4BASE)/digits+hits/hits/include \
            -I$(G4BASE)/digits+hits/digits/include \
            -I$(G4BASE)/event/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/intercoms/include

include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib

