# $Id: GNUmakefile,v 1.2 1999-04-16 09:06:07 mora Exp $
# -------------------------------------------------------------
# GNUmakefile for visualization library.  John Allison, 5/7/95.

name := G4vis_management

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4VIS_BUILD.gmk
G4TMPDIR = $(G4TMP)/$(G4SYSTEM)/$(name)

CPPFLAGS += -I$(G4BASE)/visualization/management/include
CPPFLAGS += -I$(G4BASE)/visualization/modeling/include
CPPFLAGS += -I$(G4BASE)/visualization/FukuiRenderer/include
CPPFLAGS += -I$(G4BASE)/visualization/G4Ray/include
CPPFLAGS += -I$(G4BASE)/visualization/OPACS/include
CPPFLAGS += -I$(G4BASE)/visualization/OpenGL/include
CPPFLAGS += -I$(G4BASE)/visualization/OpenInventor/include
CPPFLAGS += -I$(G4BASE)/visualization/VRML/include
CPPFLAGS += -I$(G4BASE)/graphics_reps/include
CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/global/HEPGeometry/include
CPPFLAGS += -I$(G4BASE)/global/HEPRandom/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/CSG/include
CPPFLAGS += -I$(G4BASE)/geometry/management/include
CPPFLAGS += -I$(G4BASE)/geometry/volumes/include
CPPFLAGS += -I$(G4BASE)/geometry/magneticfield/include
CPPFLAGS += -I$(G4BASE)/intercoms/include
CPPFLAGS += -I$(G4BASE)/event/include
CPPFLAGS += -I$(G4BASE)/tracking/include
CPPFLAGS += -I$(G4BASE)/track/include
CPPFLAGS += -I$(G4BASE)/digits+hits/hits/include
CPPFLAGS += -I$(G4BASE)/digits+hits/digits/include
CPPFLAGS += -I$(G4BASE)/digits+hits/detector/include
CPPFLAGS += -I$(G4BASE)/processes/management/include
CPPFLAGS += -I$(G4BASE)/particles/management/include
CPPFLAGS += -I$(G4BASE)/materials/include

include $(G4INSTALL)/config/common.gmk

ifeq ($(G4SYSTEM),WIN32-VC-NICE)
CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS))
$(G4TMPDIR)/G4VisManMessSet.o: src/G4VisManMessSet.cc
	$(CXX) $(CXXFLAGS_WITHOUT_O) $(CPPFLAGS) -c $(OUT_OBJ)$@ src/G4VisManMessSet.cc
endif
