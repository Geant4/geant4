# $Id: GNUmakefile 66892 2013-01-17 10:57:59Z gunter $
# ------------------------------------------------------------
# GNUmakefile for scoring process library.  
# ------------------------------------------------------------

name := G4scoring

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += \
	    -I$(G4BASE)/global/management/include \
	    -I$(G4BASE)/global/HEPRandom/include \
	    -I$(G4BASE)/global/HEPGeometry/include \
	    -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/navigation/include \
	    -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/geometry/biasing/include \
	    -I$(G4BASE)/track/include \
	    -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/electromagnetic/utils/include \
            -I$(G4BASE)/processes/electromagnetic/muons/include \
            -I$(G4BASE)/processes/cuts/include \
	    -I$(G4BASE)/particles/management/include \
	    -I$(G4BASE)/materials/include \
	    -I$(G4BASE)/intercoms/include \
	    -I$(G4BASE)/digits_hits/detector/include \
	    -I$(G4BASE)/digits_hits/hits/include \
	    -I$(G4BASE)/processes/transportation/include

include $(G4INSTALL)/config/common.gmk
