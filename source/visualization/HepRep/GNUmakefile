# $Id: GNUmakefile,v 1.14 2003-11-17 19:32:34 duns Exp $

name := G4visHepRep

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4VIS_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

CPPFLAGS += -DDEBUG -DZ_PREFIX
CPPFLAGS += -I$(G4BASE)/visualization/HepRep/include/ref
CPPFLAGS += -I$(G4BASE)/visualization/HepRep/include/zipios++
CPPFLAGS += -I$(G4BASE)/visualization/HepRep/include/zlib

ifeq ($(__GNUC__),2)
ifeq ($(__GNUC_MINOR__),95)
ifeq ($(__GNUC_PATCHLEVEL__),3)
# Add extra directory for sstream 2.95.2 (is available in 2.95.3)
CPPFLAGS += -I$(G4BASE)/visualization/HepRep/include/g++iostream
endif
endif
endif

CPPFLAGS += -I$(G4BASE)/digits_hits/hits/include
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
