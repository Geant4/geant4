# $Id: GNUmakefile,v 1.1 1999-01-07 16:08:00 gunter Exp $
# -------------------------------------------------------------
# GNUmakefile for STEP library.  Gabriele Cosmo, 15/11/96.
# -------------------------------------------------------------

name := G4step

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS  += \
             -I$(G4BASE)/global/management/include \
             -I$(G4BASE)/global/HEPRandom/include \
             -I$(G4BASE)/global/HEPGeometry/include

include $(G4INSTALL)/config/common.gmk

# Fix for ISO AP203 geometry standard files to be always compiled without optimization.

CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))

$(G4TMP)/$(G4SYSTEM)/$(name)/SdaiCONFIG_CONTROL_DESIGN.o: src/SdaiCONFIG_CONTROL_DESIGN.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/SdaiCONFIG_CONTROL_DESIGN.cc
$(G4TMP)/$(G4SYSTEM)/$(name)/SdaiCONFIG_CONTROL_DESIGN_init.o: src/SdaiCONFIG_CONTROL_DESIGN_init.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/SdaiCONFIG_CONTROL_DESIGN_init.cc
