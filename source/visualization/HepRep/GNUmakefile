# $Id: GNUmakefile,v 1.8 2003-07-11 16:02:36 duns Exp $

name := G4visHepRep

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4VIS_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

CPPFLAGS += -DDEBUG -DZ_PREFIX -D_GNU_SOURCE
CPPFLAGS += -I$(G4BASE)/visualization/HepRep/include/ref
CPPFLAGS += -I$(G4BASE)/visualization/HepRep/include/zipios++
CPPFLAGS += -I$(G4BASE)/visualization/HepRep/include/zlib

ifneq ($(G4_NO_STD_IOSTREAM),)
# Add extra directory for sstream 2.95.2 (is available in 2.95.3)
CPPFLAGS += -I$(G4BASE)/visualization/HepRep/include/g++iostream
endif

CPPFLAGS += -I$(G4BASE)/digits+hits/hits/include
CPPFLAGS += -I$(G4BASE)/tracking/include
CPPFLAGS += -I$(G4BASE)/visualization/management/include
CPPFLAGS += -I$(G4BASE)/visualization/modeling/include
CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/global/HEPGeometry/include
CPPFLAGS += -I$(G4BASE)/graphics_reps/include
CPPFLAGS += -I$(G4BASE)/intercoms/include
CPPFLAGS += -I$(G4BASE)/geometry/management/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/CSG/include
CPPFLAGS += -I$(G4BASE)/geometry/volumes/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/specific/include
CPPFLAGS += -I$(G4BASE)/materials/include

include $(G4INSTALL)/config/common.gmk
