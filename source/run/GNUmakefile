# $Id: GNUmakefile,v 1.6 2002-07-11 16:18:35 dressel Exp $
# ----------------------------------------------------------
# GNUmakefile for run library.  Gabriele Cosmo, 3/4/1997.
# ----------------------------------------------------------

name := G4run

ifndef G4INSTALL
  G4INSTALL = ../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/electromagnetic/utils/include \
            -I$(G4BASE)/processes/transportation/include \
            -I$(G4BASE)/processes/decay/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/tracking/include \
            -I$(G4BASE)/digits+hits/detector/include \
            -I$(G4BASE)/digits+hits/hits/include \
            -I$(G4BASE)/digits+hits/digits/include \
            -I$(G4BASE)/event/include \
            -I$(G4BASE)/intercoms/include \
	    -I$(G4BASE)/geometry/biasing/include 

include $(G4INSTALL)/config/common.gmk

.PHONY: global
 
global: lib

