# $Id: GNUmakefile,v 1.2 1999-04-15 15:28:00 johna Exp $

# lib name
name := G4VRML

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4VIS_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

CPPFLAGS += -I$(G4BASE)/visualization/management/include
CPPFLAGS += -I$(G4BASE)/visualization/modeling/include
CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/global/HEPGeometry/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/CSG/include
CPPFLAGS += -I$(G4BASE)/graphics_reps/include
CPPFLAGS += -I$(G4BASE)/intercoms/include
CPPFLAGS += -I$(G4BASE)/geometry/management/include
CPPFLAGS += -I$(G4BASE)/geometry/volumes/include
CPPFLAGS += -I$(G4BASE)/visualization/VRML/include

ifdef G4VIS_BUILD_VRML_DRIVER
  CPPFLAGS += -DVRML
endif

ifdef G4VIS_BUILD_VRMLFILE_DRIVER
  CPPFLAGS += -DVRMLFILE
endif

include $(G4INSTALL)/config/common.gmk
