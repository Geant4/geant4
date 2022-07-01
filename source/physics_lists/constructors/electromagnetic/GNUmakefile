# ---------------------------------------------------------------------------
# GNUmakefile for physics_lists/constructors/electromagnetic library.  
#  Gunter Folger 10-Jan-2012.
# ---------------------------------------------------------------------------

name := G4phys_ctor_em

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk


CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPNumerics/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/navigation/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/particles/shortlived/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/run/include \
	    -I$(G4BASE)/physics_lists/builders/include \
	    -I$(G4BASE)/physics_lists/constructors/factory/include \
	    -I$(G4BASE)/physics_lists/util/include \
	    -I$(G4BASE)/processes/hadronic/cross_sections/include \
	    -I$(G4BASE)/processes/hadronic/management/include \
	    -I$(G4BASE)/processes/hadronic/processes/include \
	    -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/transportation/include \
            -I$(G4BASE)/processes/decay/include \
            -I$(G4BASE)/processes/cuts/include \
            -I$(G4BASE)/processes/optical/include \
            -I$(G4BASE)/processes/electromagnetic/highenergy/include \
            -I$(G4BASE)/processes/electromagnetic/standard/include \
            -I$(G4BASE)/processes/electromagnetic/muons/include \
            -I$(G4BASE)/processes/electromagnetic/utils/include \
            -I$(G4BASE)/processes/electromagnetic/xrays/include \
            -I$(G4BASE)/processes/electromagnetic/lowenergy/include \
            -I$(G4BASE)/processes/electromagnetic/dna/management/include \
            -I$(G4BASE)/processes/electromagnetic/dna/molecules/management/include \
            -I$(G4BASE)/processes/electromagnetic/dna/molecules/types/include \
            -I$(G4BASE)/processes/electromagnetic/dna/models/include \
            -I$(G4BASE)/processes/electromagnetic/dna/processes/include \
            -I$(G4BASE)/processes/electromagnetic/dna/utils/include \
            -I$(G4BASE)/materials/include

include $(G4INSTALL)/config/common.gmk
