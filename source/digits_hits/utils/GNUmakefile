# $Id: GNUmakefile 99169 2016-09-07 09:29:14Z gcosmo $
# --------------------------------------------------------------
# GNUmakefile for digits+hits library.  Makoto Asai, 1/11/96.
# --------------------------------------------------------------

name := G4detutils

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4DIGI_ALLOC_EXPORT
CPPFLAGS += \
	    -I$(G4BASE)/global/HEPRandom/include \
	    -I$(G4BASE)/global/management/include \
	    -I$(G4BASE)/global/HEPGeometry/include \
	    -I$(G4BASE)/global/HEPNumerics/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/materials/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/divisions/include \
            -I$(G4BASE)/geometry/navigation/include \
	    -I$(G4BASE)/geometry/solids/CSG/include \
            -I$(G4BASE)/digits_hits/hits/include \
            -I$(G4BASE)/digits_hits/detector/include \
            -I$(G4BASE)/digits_hits/digits/include \
            -I$(G4BASE)/digits_hits/scorer/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/graphics_reps/include

include $(G4INSTALL)/config/common.gmk

