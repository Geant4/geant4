# $Id: GNUmakefile 81945 2014-06-06 16:00:38Z gcosmo $
# ----------------------------------------------------------------
# GNUmakefile for hadronic management library.  G.Folger 10-Dec-97
# ----------------------------------------------------------------


name := G4hadronic_util

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

ifdef PRECOMPOUND_TEST
  CPPFLAGS += -DPRECOMPOUND_TEST
endif

CPPFLAGS += -DG4HADRONIC_ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/materials/include

include $(G4INSTALL)/config/common.gmk

ifneq (,$(findstring WIN32-VC,$(G4SYSTEM)))

  CXXFLAGS_WITH_ZM = $(filter-out -Zm%,$(CXXFLAGS)) -Zm500

$(G4TMPDIR)/G4HadronicWhiteBoard.o: src/G4HadronicWhiteBoard.cc
	@$(ECHO) Compiling G4HadronicWhiteBoard.cc ...
	@$(CXX) $(CXXFLAGS_WITH_ZM) $(CPPFLAGS) -c $(OUT_OBJ)$@ \
	  src/G4HadronicWhiteBoard.cc

else

  CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
  CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))

$(G4TMPDIR)/G4HadronicWhiteBoard.o: src/G4HadronicWhiteBoard.cc
	@$(ECHO) Compiling G4HadronicWhiteBoard.cc ...
	@$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c -o $@ \
	  src/G4HadronicWhiteBoard.cc


endif
