
name := G4ToolsSG

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
include $(G4INSTALL)/config/interactivity.gmk

CPPFLAGS += -I$(G4BASE)/visualization/management/include
CPPFLAGS += -I$(G4BASE)/visualization/modeling/include
CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/global/HEPGeometry/include
CPPFLAGS += -I$(G4BASE)/graphics_reps/include
CPPFLAGS += -I$(G4BASE)/intercoms/include
CPPFLAGS += -I$(G4BASE)/geometry/management/include
CPPFLAGS += -I$(G4BASE)/geometry/navigation/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/CSG/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/specific/include
CPPFLAGS += -I$(G4BASE)/geometry/volumes/include
CPPFLAGS += -I$(G4BASE)/tracking/include
CPPFLAGS += -I$(G4BASE)/digits_hits/hits/include

# Locally adjust source list before including common.gmk
sources := src/G4ToolsSGSceneHandler.cc

CPPFLAGS += -I$(G4BASE)/externals/g4tools/include
sources += src/G4ToolsSGOffscreen.cc

ifdef G4VIS_BUILD_TOOLSSG_X11_GLES_DRIVER
  sources += src/G4ToolsSGX11GLES.cc
endif

ifdef G4VIS_BUILD_TOOLSSG_WINDOWS_GLES_DRIVER
  sources += src/G4ToolsSGWindowsGLES.cc
endif

ifdef G4VIS_BUILD_TOOLSSG_XT_GLES_DRIVER
  sources += src/G4ToolsSGXtGLES.cc
endif

ifdef G4VIS_BUILD_TOOLSSG_QT_GLES_DRIVER
  sources += src/G4ToolsSGQtGLES.cc src/G4ToolsSGQtViewer.cc
endif

include $(G4INSTALL)/config/moc.gmk
include $(G4INSTALL)/config/common.gmk
