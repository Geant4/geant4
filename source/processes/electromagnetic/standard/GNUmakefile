# $Id: GNUmakefile,v 1.3 2000-01-07 13:35:06 grichine Exp $
# --------------------------------------------------------------------
# GNUmakefile for electromagnetic sub-library.  John Allison, 25/6/98.
# --------------------------------------------------------------------

name := G4emstandard

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
            -I$(G4BASE)/processes/electromagnetic/muons/include \
            -I$(G4BASE)/management/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/materials/include

include $(G4INSTALL)/config/common.gmk

CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))

ifeq ($(G4SYSTEM),HP-aCC)
COMPILER := $(shell aCC -V 2>&1)
ifeq ($(COMPILER), aCC: HP ANSI C++ B3910B A.01.15)
$(G4TMP)/$(G4SYSTEM)/$(name)/G4IeEnergyLoss.o: src/G4IeEnergyLoss.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/G4IeEnergyLoss.cc
$(G4TMP)/$(G4SYSTEM)/$(name)/G4IhEnergyLoss.o: src/G4IhEnergyLoss.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/G4IhEnergyLoss.cc
$(G4TMP)/$(G4SYSTEM)/$(name)/G4eEnergyLoss.o: src/G4eEnergyLoss.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/G4eEnergyLoss.cc
$(G4TMP)/$(G4SYSTEM)/$(name)/G4eEnergyLossPlus.o: src/G4eEnergyLossPlus.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/G4eEnergyLossPlus.cc
$(G4TMP)/$(G4SYSTEM)/$(name)/G4hEnergyLoss.o: src/G4hEnergyLoss.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/G4hEnergyLoss.cc
$(G4TMP)/$(G4SYSTEM)/$(name)/G4hEnergyLossPlus.o: src/G4hEnergyLossPlus.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/G4hEnergyLossPlus.cc
endif
endif
