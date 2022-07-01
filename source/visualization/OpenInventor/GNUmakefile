
name := G4OpenInventor

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

GLOBLIBS  = libG4vis_management.lib libG4modeling.lib libG4tools.lib
ifdef G4LIB_BUILD_ZLIB
  GLOBLIBS += libG4zlib.lib
endif
GLOBLIBS += libG4interfaces.lib libG4run.lib libG4event.lib
GLOBLIBS += libG4tracking.lib libG4processes.lib libG4digits_hits.lib
GLOBLIBS += libG4track.lib libG4particles.lib libG4geometry.lib
GLOBLIBS += libG4materials.lib libG4graphics_reps.lib
GLOBLIBS += libG4intercoms.lib libG4global.lib

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4VIS_BUILD.gmk
include $(G4INSTALL)/config/G4UI_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

CPPFLAGS += -I$(G4BASE)/visualization/management/include
CPPFLAGS += -I$(G4BASE)/visualization/modeling/include
CPPFLAGS += -I$(G4BASE)/externals/g4tools/include
CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/global/HEPGeometry/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/CSG/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/specific/include
CPPFLAGS += -I$(G4BASE)/graphics_reps/include
CPPFLAGS += -I$(G4BASE)/geometry/management/include
CPPFLAGS += -I$(G4BASE)/interfaces/basic/include
CPPFLAGS += -I$(G4BASE)/interfaces/common/include
CPPFLAGS += -I$(G4BASE)/intercoms/include
CPPFLAGS += -I$(G4BASE)/materials/include
CPPFLAGS += -I$(G4BASE)/global/HEPRandom/include
CPPFLAGS += -I$(G4BASE)/tracking/include
CPPFLAGS += -I$(G4BASE)/digits_hits/hits/include

# Locally adjust source list before including common.gmk
sources := $(addprefix src/, \
  G4OpenInventor.cc \
  G4OpenInventorSceneHandler.cc \
  G4OpenInventorTransform3D.cc \
  G4OpenInventorViewer.cc \
  G4VisFeaturesOfOpenInventor.cc \
  SbPainter.cc \
  SbPainterPS.cc \
  SoAlternateRepAction.cc \
  SoBox.cc \
  SoCons.cc \
  SoCounterAction.cc \
  SoDetectorTreeKit.cc \
  SoGL2PSAction.cc \
  SoImageWriter.cc \
  SoMarkerSet.cc \
  SoPolyhedron.cc \
  SoStyleCache.cc \
  SoTrap.cc \
  SoTrd.cc \
  SoTubs.cc)

ifdef G4VIS_BUILD_OIX_DRIVER
  sources += $(addprefix src/, \
    G4OpenInventorXt.cc \
    G4OpenInventorXtExaminerViewer.cc \
    G4OpenInventorXtExaminerViewerMessenger.cc \
    G4OpenInventorXtExtended.cc \
    G4OpenInventorXtExtendedViewer.cc \
    G4OpenInventorXtViewer.cc \
    wheelmouse.cc)
endif

ifdef G4VIS_BUILD_OIQT_DRIVER
  sources += $(addprefix src/, \
    G4OpenInventorQt.cc \
    G4OpenInventorQtExaminerViewer.cc \
    G4OpenInventorQtViewer.cc \
    G4SoQt.cc)
endif

ifdef G4VIS_BUILD_OIWIN32_DRIVER
  sources += $(addprefix src/, \
    G4OpenInventorWin.cc \
    G4OpenInventorWinViewer.cc)
endif

include $(G4INSTALL)/config/common.gmk

debug:
	@echo $(CPPFLAGS)
