# $Id: GNUmakefile,v 1.2 2003/10/03 10:29:19 gcosmo Exp $
# --------------------------------------------------------------
# GNUmakefile for digits+hits library.  Makoto Asai, 1/11/96.
# --------------------------------------------------------------

name := G4detector

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += \
	    -I$(G4BASE)/global/HEPRandom/include \
	    -I$(G4BASE)/global/management/include \
	    -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/navigation/include \
	    -I$(G4BASE)/geometry/biasing/include \
            -I$(G4BASE)/digits_hits/hits/include \
            -I$(G4BASE)/intercoms/include

include $(G4INSTALL)/config/common.gmk

