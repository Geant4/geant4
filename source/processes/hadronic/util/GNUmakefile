# $Id: GNUmakefile,v 1.1 1999-01-07 16:13:48 gunter Exp $
# ----------------------------------------------------------------
# GNUmakefile for hadronic management library.  G.Folger 10-Dec-97
# ----------------------------------------------------------------

name := G4hadronic_util

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

G4TMPDIR = $(G4TMP)/$(G4SYSTEM)/$(name)

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

# Temporary fix for HP-CC compiler. It cannot handle full optimisation for
# those files.

ifneq (,$(findstring HP,$(G4SYSTEM)))
    CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS))

$(G4TMPDIR)/G4NucleiPropertiesTable.o: src/G4NucleiPropertiesTable.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c -o $@ src/G4NucleiPropertiesTable.cc
$(G4TMPDIR)/G4ReactionDynamics.o: src/G4ReactionDynamics.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c -o $@ src/G4ReactionDynamics.cc
endif

ifneq (,$(findstring SUN,$(G4SYSTEM)))
    CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
$(G4TMPDIR)/G4ReactionDynamics.o: src/G4ReactionDynamics.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c -o $@ src/G4ReactionDynamics.cc
endif

ifneq (,$(findstring WIN32-VC,$(G4SYSTEM)))
    CXXFLAGS_WITH_ZM = $(CXXFLAGS) -Zm500

$(G4TMPDIR)/G4NucleiPropertiesTable.o: src/G4NucleiPropertiesTable.cc
	$(CXX) $(CXXFLAGS_WITH_ZM) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/G4NucleiPropertiesTable.cc
endif
