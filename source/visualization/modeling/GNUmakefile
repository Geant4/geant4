# $Id: GNUmakefile,v 1.6 2003/11/06 15:15:08 johna Exp $
# --------------------------------------------------------------
# GNUmakefile for modeling library.  John Allison, 31/12/1997.
# --------------------------------------------------------------

name := G4modeling

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk

CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/graphics_reps/include
CPPFLAGS += -I$(G4BASE)/global/HEPGeometry/include
CPPFLAGS += -I$(G4BASE)/geometry/management/include
CPPFLAGS += -I$(G4BASE)/geometry/volumes/include
CPPFLAGS += -I$(G4BASE)/geometry/navigation/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/CSG/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/specific/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/Boolean/include
CPPFLAGS += -I$(G4BASE)/materials/include
CPPFLAGS += -I$(G4BASE)/global/HEPRandom/include
CPPFLAGS += -I$(G4BASE)/processes/management/include
CPPFLAGS += -I$(G4BASE)/particles/management/include
CPPFLAGS += -I$(G4BASE)/run/include
CPPFLAGS += -I$(G4BASE)/event/include
CPPFLAGS += -I$(G4BASE)/digits_hits/hits/include
CPPFLAGS += -I$(G4BASE)/digits_hits/digits/include
CPPFLAGS += -I$(G4BASE)/tracking/include
CPPFLAGS += -I$(G4BASE)/track/include
CPPFLAGS += -I$(G4BASE)/intercoms/include

include $(G4INSTALL)/config/common.gmk
