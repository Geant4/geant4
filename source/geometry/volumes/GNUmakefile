# $Id: GNUmakefile,v 1.4 2001-10-24 15:33:47 gcosmo Exp $
# --------------------------------------------------------------------
# GNUmakefile for geometry/volumes library.  Gabriele Cosmo, 16/11/96.
# --------------------------------------------------------------------

name := G4volumes

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/graphics_reps/include \
            -I$(G4BASE)/intercoms/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
            -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/geometry/verification/include \
            -I$(G4BASE)/geometry/management/include

include $(G4INSTALL)/config/common.gmk

CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))

ifeq ($(G4SYSTEM),SUN-CC4)
$(G4TMP)/$(G4SYSTEM)/$(name)/G4VoxelNavigation.o: src/G4VoxelNavigation.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/G4VoxelNavigation.cc
endif

