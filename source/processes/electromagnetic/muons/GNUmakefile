# $Id: GNUmakefile,v 1.3 1999-06-05 13:20:23 stesting Exp $
# --------------------------------------------------------------------
# GNUmakefile for electromagnetic sub-library.  John Allison, 25/6/98.
# --------------------------------------------------------------------

name := G4muons

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/electromagnetic/utils/include \
            -I$(G4BASE)/processes/electromagnetic/standard/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/materials/include

include $(G4INSTALL)/config/common.gmk

CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))

ifeq ($(G4SYSTEM),HP-aCC)
COMPILER := $(shell aCC -V 2>&1)
ifeq ($(COMPILER), aCC: HP ANSI C++ B3910B A.01.15)
$(G4TMP)/$(G4SYSTEM)/$(name)/G4MuEnergyLoss.o: src/G4MuEnergyLoss.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/G4MuEnergyLoss.cc
endif
endif
