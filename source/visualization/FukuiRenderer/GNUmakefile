# $Id: GNUmakefile,v 1.1 1999-01-07 16:14:34 gunter Exp $

name := G4FR

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4VIS_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

CPPFLAGS += -I$(G4BASE)/visualization/management/include
CPPFLAGS += -I$(G4BASE)/visualization/modeling/include
CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/graphics_reps/include
CPPFLAGS += -I$(G4BASE)/global/HEPGeometry/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/CSG/include
CPPFLAGS += -I$(G4BASE)/geometry/management/include
CPPFLAGS += -I$(G4BASE)/geometry/volumes/include

include $(G4INSTALL)/config/common.gmk
