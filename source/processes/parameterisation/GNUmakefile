# $Id: GNUmakefile,v 1.4 2000-07-14 08:18:04 grichine Exp $
# ------------------------------------------------------------
# GNUmakefile for transportation library.  G.Folger 10-Dec-97.
# ------------------------------------------------------------

name := G4parameterisation

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
	    -I$(G4BASE)/global/HEPRandom/include \
	    -I$(G4BASE)/global/HEPGeometry/include \
	    -I$(G4BASE)/global/HEPNumerics/include \
	    -I$(G4BASE)/geometry/management/include \
	    -I$(G4BASE)/geometry/volumes/include \
	    -I$(G4BASE)/geometry/magneticfield/include \
	    -I$(G4BASE)/geometry/solids/CSG/include \
	    -I$(G4BASE)/intercoms/include   \
	    -I$(G4BASE)/track/include \
	    -I$(G4BASE)/tracking/include \
	    -I$(G4BASE)/step/include \
	    -I$(G4BASE)/digits+hits/detector/include \
	    -I$(G4BASE)/digits+hits/digits/include \
	    -I$(G4BASE)/digits+hits/hits/include \
	    -I$(G4BASE)/processes/management/include \
	    -I$(G4BASE)/processes/electromagnetic/standard/include \
	    -I$(G4BASE)/processes/electromagnetic/utils/include \
	    -I$(G4BASE)/particles/management/include \
	    -I$(G4BASE)/particles/bosons/include \
	    -I$(G4BASE)/particles/hadrons/include \
	    -I$(G4BASE)/particles/hadrons/barions/include \
	    -I$(G4BASE)/particles/leptons/include \
	    -I$(G4BASE)/materials/include \

include $(G4INSTALL)/config/common.gmk
