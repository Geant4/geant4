# -----------------------------------------------------------
# GNUmakefile for hadronic library.  JPW, 2003.
# -----------------------------------------------------------

name := G4had_string_diff

ifndef G4INSTALL
  G4INSTALL = ../../../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4HADRONIC_ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPNumerics/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/hadronic/management/include/ \
            -I$(G4BASE)/processes/hadronic/cross_sections/include/ \
            -I$(G4BASE)/processes/hadronic/util/include \
	    -I$(G4BASE)/processes/hadronic/models/parton_string/diffraction/include \
	    -I$(G4BASE)/processes/hadronic/models/parton_string/management/include \
	    -I$(G4BASE)/processes/hadronic/models/parton_string/hadronization/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/particles/shortlived/include \
            -I$(G4BASE)/materials/include
           
include $(G4INSTALL)/config/common.gmk

