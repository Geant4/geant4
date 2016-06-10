# $Id: GNUmakefile 66892 2013-01-17 10:57:59Z gunter $
# ----------------------------------------------------------------
# GNUmakefile for hadronic management library.  G.Folger 10-Dec-97
# ----------------------------------------------------------------

name := G4hadronic_mgt

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -DG4HADRONIC_ALLOC_EXPORT
CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/navigation/include \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/processes/hadronic/models/management/include \
            -I$(G4BASE)/processes/hadronic/cross_sections/include \
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

else 

  CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
  CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))

$(G4TMPDIR)/G4HadronicInteractionWrapper.o: src/G4HadronicInteractionWrapper.cc
	@$(ECHO) Compiling G4HadronicInteractionWrapper.cc ...
	@$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -g -c -o $@ \
	  src/G4HadronicInteractionWrapper.cc
endif
