# $Id: GNUmakefile,v 1.1 1999-01-07 16:12:25 gunter Exp $
# -----------------------------------------------------------
# GNUmakefile for hadronic library.  Gabriele Cosmo, 18/9/96.
# -----------------------------------------------------------

name := G4hadronic_HE

ifndef G4INSTALL
  G4INSTALL = ../../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

G4TMPDIR = $(G4TMP)/$(G4SYSTEM)/$(name)

CPPFLAGS += -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/track/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/processes/management/include \
            -I$(G4BASE)/processes/hadronic/management/include/ \
            -I$(G4BASE)/processes/hadronic/util/include \
            -I$(G4BASE)/particles/management/include \
            -I$(G4BASE)/particles/leptons/include \
            -I$(G4BASE)/particles/bosons/include \
            -I$(G4BASE)/particles/hadrons/mesons/include \
            -I$(G4BASE)/particles/hadrons/barions/include \
            -I$(G4BASE)/particles/hadrons/ions/include \
            -I$(G4BASE)/materials/include
           
include $(G4INSTALL)/config/common.gmk

# Temporary fix for HP-CC compiler. It cannot handle full optimisation for
# the specified files ...

ifneq (,$(findstring HP,$(G4SYSTEM)))
    CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS))

$(G4TMPDIR)/G4GHEInelastic.o: src/G4GHEInelastic.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c -o $@ src/G4GHEInelastic.cc

$(G4TMPDIR)/G4GHEInelasticModel.o: src/G4GHEInelasticModel.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c -o $@ src/G4GHEInelasticModel.cc
endif
