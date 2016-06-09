# $Id: GNUmakefile,v 1.23 2005/06/28 08:47:50 gcosmo Exp $

name := G4visHepRep

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

GLOBLIBS  = libG4vis_management.lib libG4modeling.lib
ifdef G4LIB_BUILD_ZLIB
  GLOBLIBS += libG4zlib.lib
endif
GLOBLIBS += libG4run.lib libG4event.lib
GLOBLIBS += libG4tracking.lib libG4processes.lib libG4digits_hits.lib
GLOBLIBS += libG4track.lib libG4particles.lib libG4geometry.lib
GLOBLIBS += libG4materials.lib libG4graphics_reps.lib
GLOBLIBS += libG4intercoms.lib libG4global.lib

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4VIS_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

ifdef G4LIB_USE_ZLIB
  CPPFLAGS += -DG4LIB_USE_ZLIB
  ifdef G4LIB_BUILD_ZLIB
    CPPFLAGS += -I$(G4BASE)/visualization/externals/zlib/include
  endif
else
  CPPFLAGS += -DCHEPREP_NO_ZLIB
endif

CPPFLAGS += -I$(G4BASE)/event/include
CPPFLAGS += -I$(G4BASE)/track/include
CPPFLAGS += -I$(G4BASE)/run/include
CPPFLAGS += -I$(G4BASE)/digits_hits/hits/include
CPPFLAGS += -I$(G4BASE)/particles/management/include
CPPFLAGS += -I$(G4BASE)/tracking/include
CPPFLAGS += -I$(G4BASE)/visualization/management/include
CPPFLAGS += -I$(G4BASE)/visualization/modeling/include
CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/global/HEPRandom/include
CPPFLAGS += -I$(G4BASE)/global/HEPGeometry/include
CPPFLAGS += -I$(G4BASE)/graphics_reps/include
CPPFLAGS += -I$(G4BASE)/intercoms/include
CPPFLAGS += -I$(G4BASE)/processes/management/include
CPPFLAGS += -I$(G4BASE)/geometry/management/include
CPPFLAGS += -I$(G4BASE)/geometry/navigation/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/CSG/include
CPPFLAGS += -I$(G4BASE)/geometry/volumes/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/specific/include
CPPFLAGS += -I$(G4BASE)/materials/include

include $(G4INSTALL)/config/common.gmk
