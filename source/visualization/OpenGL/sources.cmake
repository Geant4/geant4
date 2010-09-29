#------------------------------------------------------------------------------
# sources.cmake
# Module : G4OpenGL
# Package: Geant4.src.G4visualization.G4OpenGL
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake,v 1.1 2010-09-29 19:12:14 bmorgan Exp $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/specific/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/interfaces/basic/include)
include_directories(${CMAKE_SOURCE_DIR}/source/interfaces/common/include)
include_directories(${CMAKE_SOURCE_DIR}/source/tracking/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/externals/gl2ps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/modeling/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4OpenGL
    HEADERS
        G4OpenGL.hh
        G4OpenGLBitMapStore.hh
        G4OpenGLFontBaseStore.hh
        G4OpenGLImmediateQt.hh
        G4OpenGLImmediateQtViewer.hh
        G4OpenGLImmediateSceneHandler.hh
        G4OpenGLImmediateViewer.hh
        G4OpenGLImmediateWin32.hh
        G4OpenGLImmediateWin32Viewer.hh
        G4OpenGLImmediateX.hh
        G4OpenGLImmediateXViewer.hh
        G4OpenGLImmediateXm.hh
        G4OpenGLImmediateXmViewer.hh
        G4OpenGLQtExportDialog.hh
        G4OpenGLQtMovieDialog.hh
        G4OpenGLQtViewer.hh
        G4OpenGLSceneHandler.hh
        G4OpenGLSceneHandler.icc
        G4OpenGLStoredQt.hh
        G4OpenGLStoredQtViewer.hh
        G4OpenGLStoredSceneHandler.hh
        G4OpenGLStoredViewer.hh
        G4OpenGLStoredWin32.hh
        G4OpenGLStoredWin32Viewer.hh
        G4OpenGLStoredX.hh
        G4OpenGLStoredXViewer.hh
        G4OpenGLStoredXm.hh
        G4OpenGLStoredXmViewer.hh
        G4OpenGLTransform3D.hh
        G4OpenGLViewer.hh
        G4OpenGLViewerMessenger.hh
        G4OpenGLWin32Viewer.hh
        G4OpenGLXViewer.hh
        G4OpenGLXmBox.hh
        G4OpenGLXmFourArrowButtons.hh
        G4OpenGLXmFramedBox.hh
        G4OpenGLXmPushButton.hh
        G4OpenGLXmRadioButton.hh
        G4OpenGLXmResources.hh
        G4OpenGLXmSeparator.hh
        G4OpenGLXmSliderBar.hh
        G4OpenGLXmTextField.hh
        G4OpenGLXmTopLevelShell.hh
        G4OpenGLXmVWidgetComponent.hh
        G4OpenGLXmVWidgetContainer.hh
        G4OpenGLXmVWidgetObject.hh
        G4OpenGLXmVWidgetShell.hh
        G4OpenGLXmViewer.hh
        G4OpenGLXmViewerMessenger.hh
    SOURCES
        G4OpenGLBitMapStore.cc
        G4OpenGLFontBaseStore.cc
        G4OpenGLImmediateQt.cc
        G4OpenGLImmediateQtViewer.cc
        G4OpenGLImmediateSceneHandler.cc
        G4OpenGLImmediateViewer.cc
        G4OpenGLImmediateWin32.cc
        G4OpenGLImmediateWin32Viewer.cc
        G4OpenGLImmediateX.cc
        G4OpenGLImmediateXViewer.cc
        G4OpenGLImmediateXm.cc
        G4OpenGLImmediateXmViewer.cc
        G4OpenGLQtExportDialog.cc
        G4OpenGLQtMovieDialog.cc
        G4OpenGLQtViewer.cc
        G4OpenGLSceneHandler.cc
        G4OpenGLStoredQt.cc
        G4OpenGLStoredQtViewer.cc
        G4OpenGLStoredSceneHandler.cc
        G4OpenGLStoredViewer.cc
        G4OpenGLStoredWin32.cc
        G4OpenGLStoredWin32Viewer.cc
        G4OpenGLStoredX.cc
        G4OpenGLStoredXViewer.cc
        G4OpenGLStoredXm.cc
        G4OpenGLStoredXmViewer.cc
        G4OpenGLTransform3D.cc
        G4OpenGLViewer.cc
        G4OpenGLViewerMessenger.cc
        G4OpenGLWin32Viewer.cc
        G4OpenGLXViewer.cc
        G4OpenGLXmBox.cc
        G4OpenGLXmConvenienceRoutines.cc
        G4OpenGLXmFourArrowButtons.cc
        G4OpenGLXmFramedBox.cc
        G4OpenGLXmMainMenubarCallbacks.cc
        G4OpenGLXmPanningCallbacks.cc
        G4OpenGLXmPushButton.cc
        G4OpenGLXmRadioButton.cc
        G4OpenGLXmRotationCallbacks.cc
        G4OpenGLXmSeparator.cc
        G4OpenGLXmSliderBar.cc
        G4OpenGLXmStyleCallbacks.cc
        G4OpenGLXmTextField.cc
        G4OpenGLXmTopLevelShell.cc
        G4OpenGLXmVWidgetComponent.cc
        G4OpenGLXmVWidgetContainer.cc
        G4OpenGLXmVWidgetObject.cc
        G4OpenGLXmVWidgetShell.cc
        G4OpenGLXmViewer.cc
        G4OpenGLXmViewerMessenger.cc
        G4OpenGLXmWindowHandlingCallbacks.cc
    GRANULAR_DEPENDENCIES
        G4UIbasic
        G4UIcommon
        G4csg
        G4geometrymng
        G4gl2ps
        G4globman
        G4graphics_reps
        G4hits
        G4intercoms
        G4modeling
        G4specsolids
        G4tracking
        G4vis_management
    GLOBAL_DEPENDENCIES
        G4digits_hits
        G4geometry
        G4gl2ps
        G4global
        G4graphics_reps
        G4intercoms
        G4interfaces
        G4modeling
        G4tracking
        G4vis_management
    LINK_LIBRARIES
)

# List any source specific properties here

