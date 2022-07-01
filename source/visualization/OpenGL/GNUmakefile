
name := G4OpenGL

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

# Definition of macro for moc files (for Qt use only)
ifdef G4VIS_BUILD_OPENGLQT_DRIVER
MOC_MACRO = -DG4VIS_BUILD_OPENGLQT_DRIVER
endif

# For debug mode
# CPPFLAGS += -DG4DEBUG_VIS_OGL

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

CPPFLAGS += -I$(G4BASE)/visualization/externals/gl2ps/include
CPPFLAGS += -I$(G4BASE)/visualization/management/include
CPPFLAGS += -I$(G4BASE)/visualization/modeling/include
CPPFLAGS += -I$(G4BASE)/externals/g4tools/include
CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/global/HEPGeometry/include
CPPFLAGS += -I$(G4BASE)/global/HEPRandom/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/CSG/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/specific/include
CPPFLAGS += -I$(G4BASE)/graphics_reps/include
CPPFLAGS += -I$(G4BASE)/geometry/management/include
CPPFLAGS += -I$(G4BASE)/geometry/navigation/include
CPPFLAGS += -I$(G4BASE)/geometry/volumes/include
CPPFLAGS += -I$(G4BASE)/processes/management/include
CPPFLAGS += -I$(G4BASE)/intercoms/include
CPPFLAGS += -I$(G4BASE)/interfaces/common/include
CPPFLAGS += -I$(G4BASE)/interfaces/basic/include
CPPFLAGS += -I$(G4BASE)/materials/include
CPPFLAGS += -I$(G4BASE)/track/include
CPPFLAGS += -I$(G4BASE)/tracking/include
CPPFLAGS += -I$(G4BASE)/digits_hits/digits/include
CPPFLAGS += -I$(G4BASE)/digits_hits/hits/include
CPPFLAGS += -I$(G4BASE)/particles/management/include
CPPFLAGS += -I$(G4BASE)/event/include
CPPFLAGS += -I$(G4BASE)/run/include
CPPFLAGS += -I$(G4BASE)/tasking/include

# Locally adjust source list before including common.gmk
sources := $(addprefix src/,G4OpenGLImmediateViewer.cc \
  G4OpenGLImmediateSceneHandler.cc \
  G4OpenGLViewer.cc \
  G4OpenGLStoredViewer.cc \
  G4OpenGLStoredSceneHandler.cc \
  G4OpenGLFontBaseStore.cc \
  G4OpenGLSceneHandler.cc \
  G4OpenGLViewerMessenger.cc \
  G4OpenGLTransform3D.cc \
  G4VisFeaturesOfOpenGL.cc \
  G4gl2ps.cc)

ifdef G4VIS_BUILD_OPENGLX_DRIVER
  sources += $(addprefix src/,\
    G4OpenGLImmediateX.cc \
    G4OpenGLImmediateXViewer.cc \
    G4OpenGLStoredX.cc \
    G4OpenGLStoredXViewer.cc \
    G4OpenGLXViewer.cc)
endif

ifdef G4VIS_BUILD_OPENGLXM_DRIVER
  sources += $(addprefix src/,\
    G4OpenGLImmediateXm.cc \
    G4OpenGLImmediateXmViewer.cc \
    G4OpenGLStoredXm.cc \
    G4OpenGLStoredXmViewer.cc \
    G4OpenGLXm.cc \
    G4OpenGLXmBox.cc \
    G4OpenGLXmConvenienceRoutines.cc \
    G4OpenGLXmFourArrowButtons.cc \
    G4OpenGLXmFramedBox.cc \
    G4OpenGLXmMainMenubarCallbacks.cc \
    G4OpenGLXmPanningCallbacks.cc \
    G4OpenGLXmPushButton.cc \
    G4OpenGLXmRadioButton.cc \
    G4OpenGLXmRotationCallbacks.cc \
    G4OpenGLXmSeparator.cc \
    G4OpenGLXmSliderBar.cc \
    G4OpenGLXmStyleCallbacks.cc \
    G4OpenGLXmTextField.cc \
    G4OpenGLXmTopLevelShell.cc \
    G4OpenGLXmViewer.cc \
    G4OpenGLXmViewerMessenger.cc \
    G4OpenGLXmVWidgetComponent.cc \
    G4OpenGLXmVWidgetContainer.cc \
    G4OpenGLXmVWidgetObject.cc \
    G4OpenGLXmVWidgetShell.cc \
    G4OpenGLXmWindowHandlingCallbacks.cc)
endif

ifdef G4VIS_BUILD_OPENGLQT_DRIVER
  sources += $(addprefix src/, \
    G4OpenGLImmediateQt.cc \
    G4OpenGLImmediateQtViewer.cc \
    G4OpenGLQt.cc \
    G4OpenGLQtExportDialog.cc \
    G4OpenGLQtMovieDialog.cc \
    G4OpenGLVboDrawer.cc \
    G4OpenGLQtViewer.cc \
    G4OpenGLStoredQt.cc \
    G4OpenGLStoredQtSceneHandler.cc \
    G4OpenGLStoredQtViewer.cc)
endif

ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER
  sources += $(addprefix src/, \
    G4OpenGLImmediateWin32.cc \
    G4OpenGLImmediateWin32Viewer.cc \
    G4OpenGLStoredWin32.cc \
    G4OpenGLStoredWin32Viewer.cc \
    G4OpenGLWin32Viewer.cc)
endif

include $(G4INSTALL)/config/moc.gmk
include $(G4INSTALL)/config/common.gmk
