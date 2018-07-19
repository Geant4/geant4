# $Id: GNUmakefile 99176 2016-09-07 09:46:36Z gcosmo $

name := G4GMocren

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

GLOBLIBS  = libG4vis_management.lib libG4modeling.lib
GLOBLIBS += libG4run.lib libG4event.lib
GLOBLIBS += libG4tracking.lib libG4processes.lib libG4digits_hits.lib
GLOBLIBS += libG4track.lib libG4particles.lib libG4geometry.lib
GLOBLIBS += libG4materials.lib libG4graphics_reps.lib
GLOBLIBS += libG4intercoms.lib libG4global.lib

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4VIS_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

# NO need OGL, then reset OGLGLAGS
OGLFLAGS =
OGLLIBS  =
# NO need QT, then reset QTGLAGS
QTFLAGS =
QTLIBS  =
GLQTLIBS =

CPPFLAGS += -I$(G4BASE)/visualization/management/include
CPPFLAGS += -I$(G4BASE)/visualization/modeling/include
CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/global/HEPGeometry/include
CPPFLAGS += -I$(G4BASE)/global/HEPNumerics/include
CPPFLAGS += -I$(G4BASE)/graphics_reps/include
CPPFLAGS += -I$(G4BASE)/intercoms/include
CPPFLAGS += -I$(G4BASE)/geometry/management/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/CSG/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/specific/include

CPPFLAGS += -I$(G4BASE)/geometry/navigation/include/

CPPFLAGS += -I$(G4BASE)/particles/management/include
CPPFLAGS += -I$(G4BASE)/event/include
CPPFLAGS += -I$(G4BASE)/tracking/include
CPPFLAGS += -I$(G4BASE)/materials/include
CPPFLAGS += -I$(G4BASE)/digits_hits/hits/include
CPPFLAGS += -I$(G4BASE)/digits_hits/digits/include
CPPFLAGS += -I$(G4BASE)/digits_hits/utils/include

#CPPFLAGS += -I$(G4BASE)/visualization/FukuiRenderer/include


#CPPFLAGS += -DG4GMocrenFileDEBUG
#PPFLAGS += -DGMOCRENIO_LITTLE_ENDIAN

ifdef G4VIS_BUILD_XXXSG_DRIVER
  CPPFLAGS += -DG4VIS_BUILD_XXXSG_DRIVER
endif

include $(G4INSTALL)/config/common.gmk
