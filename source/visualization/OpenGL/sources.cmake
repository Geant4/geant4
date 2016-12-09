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
# $Id: sources.cmake 99440 2016-09-22 08:34:04Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/digits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/event/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/specific/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/interfaces/basic/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/interfaces/common/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/run/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)
include_directories(${CMAKE_SOURCE_DIR}/source/tracking/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/externals/gl2ps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/modeling/include)

#
# Module has optional sources
#
include(Geant4MacroDefineModule)

#
# Define the core sources, includes and libraries which all Geant4
# OpenGL implementations use
#
set(G4VIS_MODULE_OPENGL_HEADERS
    G4OpenGL.hh
    G4OpenGLImmediateViewer.hh
    G4OpenGLImmediateSceneHandler.hh
    G4OpenGLViewer.hh
    G4OpenGLStoredViewer.hh
    G4OpenGLStoredSceneHandler.hh
    G4OpenGLFontBaseStore.hh
    G4OpenGLSceneHandler.hh
    G4OpenGLSceneHandler.icc
    G4OpenGLViewerMessenger.hh
    G4OpenGLTransform3D.hh
    G4VisFeaturesOfOpenGL.hh
)

set(G4VIS_MODULE_OPENGL_SOURCES
    G4OpenGLImmediateViewer.cc
    G4OpenGLImmediateSceneHandler.cc
    G4OpenGLViewer.cc
    G4OpenGLStoredViewer.cc
    G4OpenGLStoredSceneHandler.cc
    G4OpenGLFontBaseStore.cc
    G4OpenGLSceneHandler.cc
    G4OpenGLViewerMessenger.cc
    G4OpenGLTransform3D.cc
    G4VisFeaturesOfOpenGL.cc
)

#
# May need OpenGL include here
#
include_directories(${OPENGL_INCLUDE_DIR})
set(G4VIS_MODULE_OPENGL_LINK_LIBRARIES ${OPENGL_LIBRARIES})

#
# All files must have the G4VIS_BUILD_OPENGL_DRIVER definition
#
add_definitions(-DG4VIS_BUILD_OPENGL_DRIVER)


#----------------------------------------------------------------------------
# Add X11 OpenGL Support if requested
#
if(GEANT4_USE_OPENGL_X11)
    #
    # Add in the extra X11 GL sources
    #
    list(APPEND G4VIS_MODULE_OPENGL_HEADERS
        G4OpenGLImmediateX.hh
        G4OpenGLImmediateXViewer.hh
        G4OpenGLStoredX.hh
        G4OpenGLStoredXViewer.hh
        G4OpenGLXViewer.hh)

    list(APPEND G4VIS_MODULE_OPENGL_SOURCES
        G4OpenGLImmediateX.cc
        G4OpenGLImmediateXViewer.cc
        G4OpenGLStoredX.cc
        G4OpenGLStoredXViewer.cc
        G4OpenGLXViewer.cc)

    # Add the includes for X11
    include_directories(${X11_INCLUDE_DIR} ${X11_Xmu_INCLUDE_PATH})

    # Add the compile definitions needed for the X11 component
    add_definitions(-DG4VIS_BUILD_OPENGLX_DRIVER)

    # Add in X11 libraries, plus Xmu library
    list(APPEND G4VIS_MODULE_OPENGL_LINK_LIBRARIES ${X11_LIBRARIES} ${X11_Xmu_LIBRARY})
endif()


#----------------------------------------------------------------------------
# Add Xm OpenGL Support if requested
#
if(GEANT4_USE_XM)
    #
    # Add in the extra Xm sources
    #
    list(APPEND G4VIS_MODULE_OPENGL_HEADERS
        G4OpenGLImmediateXm.hh
        G4OpenGLImmediateXmViewer.hh
        G4OpenGLStoredXm.hh
        G4OpenGLStoredXmViewer.hh
        G4OpenGLXm.hh
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
        G4OpenGLXmViewer.hh
        G4OpenGLXmViewerMessenger.hh
        G4OpenGLXmVWidgetComponent.hh
        G4OpenGLXmVWidgetContainer.hh
        G4OpenGLXmVWidgetObject.hh
        G4OpenGLXmVWidgetShell.hh
    )

    list(APPEND G4VIS_MODULE_OPENGL_SOURCES
        G4OpenGLImmediateXm.cc
        G4OpenGLImmediateXmViewer.cc
        G4OpenGLStoredXm.cc
        G4OpenGLStoredXmViewer.cc
        G4OpenGLXm.cc
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
        G4OpenGLXmViewer.cc
        G4OpenGLXmViewerMessenger.cc
        G4OpenGLXmVWidgetComponent.cc
        G4OpenGLXmVWidgetContainer.cc
        G4OpenGLXmVWidgetObject.cc
        G4OpenGLXmVWidgetShell.cc
        G4OpenGLXmWindowHandlingCallbacks.cc
    )

    # Special case of building Xm without X11
    if(NOT GEANT4_USE_OPENGL_X11)
      list(APPEND G4VIS_MODULE_OPENGL_HEADERS G4OpenGLXViewer.hh)
      list(APPEND G4VIS_MODULE_OPENGL_SOURCES G4OpenGLXViewer.cc)
      add_definitions(-DG4VIS_BUILD_OPENGLX_DRIVER)
    endif()

    # Add the includes for X11, Xmu and Motif
    include_directories(
        ${X11_INCLUDE_DIR}
        ${X11_Xmu_INCLUDE_PATH}
        ${MOTIF_INCLUDE_DIR}
    )

    # Add the compile definitions needed for the Xm component, remembering
    # to add those for the UI part as well!!!
    add_definitions(-DG4VIS_BUILD_OPENGLXM_DRIVER -DG4INTY_BUILD_XT -DG4UI_BUILD_XM_SESSION)

    # Add in Xm, X11 libraries, plus Xmu library
    set(G4VIS_MODULE_OPENGL_LINK_LIBRARIES ${MOTIF_LIBRARIES} ${X11_LIBRARIES} ${X11_Xmu_LIBRARY} ${G4VIS_MODULE_OPENGL_LINK_LIBRARIES})
    #list(APPEND G4VIS_MODULE_OPENGL_LINK_LIBRARIES
    #       ${MOTIF_LIBRARIES}
    #       ${X11_LIBRARIES}
    #       ${X11_Xmu_LIBRARY}
    #)
endif()


#----------------------------------------------------------------------------
# Add Qt OpenGL support if requested
#
if(GEANT4_USE_QT)
    #
    # Add in the extra Qt GL sources
    #
    list(APPEND G4VIS_MODULE_OPENGL_HEADERS
        G4OpenGLImmediateQt.hh
        G4OpenGLImmediateQtViewer.hh
        G4OpenGLQt.hh
        G4OpenGLQtExportDialog.hh
        G4OpenGLQtMovieDialog.hh
        G4OpenGLVboDrawer.hh
        G4OpenGLQtViewer.hh
        G4OpenGLStoredQt.hh
        G4OpenGLStoredQtSceneHandler.hh
        G4OpenGLStoredQtViewer.hh)

    list(APPEND G4VIS_MODULE_OPENGL_SOURCES
        G4OpenGLImmediateQt.cc
        G4OpenGLImmediateQtViewer.cc
        G4OpenGLQt.cc
        G4OpenGLQtExportDialog.cc
        G4OpenGLQtMovieDialog.cc
        G4OpenGLVboDrawer.cc
        G4OpenGLQtViewer.cc
        G4OpenGLStoredQt.cc
        G4OpenGLStoredQtSceneHandler.cc
        G4OpenGLStoredQtViewer.cc)


    # Include the UseQt file to build the moc wrappers
    include(${QT_USE_FILE})

    # Add the moc sources - must use absolute path to the files
    QT4_WRAP_CPP(G4OPENGL_MOC_SOURCES
        ${CMAKE_SOURCE_DIR}/source/visualization/OpenGL/include/G4OpenGLQtExportDialog.hh
        ${CMAKE_SOURCE_DIR}/source/visualization/OpenGL/include/G4OpenGLQtMovieDialog.hh
        ${CMAKE_SOURCE_DIR}/source/visualization/OpenGL/include/G4OpenGLQtViewer.hh
         OPTIONS -DG4VIS_BUILD_OPENGLQT_DRIVER)

    list(APPEND G4VIS_MODULE_OPENGL_SOURCES ${G4OPENGL_MOC_SOURCES})

    # Add the definitions - these will also be used to compile the moc sources
    # Argh.. Have to remember about INTY and UI because of their use...
    add_definitions(-DG4VIS_BUILD_OPENGLQT_DRIVER -DG4INTY_BUILD_QT
        -DG4UI_BUILD_QT_SESSION)

    # Add in Qt libraries
    list(APPEND G4VIS_MODULE_OPENGL_LINK_LIBRARIES ${QT_LIBRARIES})
endif()


#----------------------------------------------------------------------------
# Add Wt OpenGL support if requested
#
if(GEANT4_USE_WT)
    #
    # Add in the extra Wt GL sources
    #
    list(APPEND G4VIS_MODULE_OPENGL_HEADERS
        G4OpenGLImmediateWt.hh
        G4OpenGLImmediateWtViewer.hh
        G4OpenGLVboDrawer.hh
        G4OpenGLWtViewer.hh)

    list(APPEND G4VIS_MODULE_OPENGL_SOURCES
        G4OpenGLImmediateWt.cc
        G4OpenGLImmediateWtViewer.cc
        G4OpenGLVboDrawer.cc
        G4OpenGLWtViewer.cc)

    # Must have Wt includes...
    include_directories(${Wt_INCLUDE_DIR})
    include_directories(${Boost_INCLUDE_DIRS})

    # Add the definitions - these will also be used to compile the moc sources
    # Argh.. Have to remember about INTY and UI because of their use...
    # Use the compile definitions to avoid "signal/slot" keyword
    # mixed with Qt and boost
    add_definitions(${WT_DEFINITIONS})
    add_definitions(-DG4VIS_BUILD_OPENGLWT_DRIVER -DG4INTY_BUILD_WT
        -DG4UI_BUILD_WT_SESSION)

    # Add in Wt libraries
    list(APPEND G4VIS_MODULE_OPENGL_LINK_LIBRARIES ${Wt_LIBRARY})
endif()


#----------------------------------------------------------------------------
# Add WIN32 support, if requested
#
if(GEANT4_USE_OPENGL_WIN32)
    #
    # Add in the extra sources
    #
    list(APPEND G4VIS_MODULE_OPENGL_HEADERS
        G4OpenGLImmediateWin32.hh
        G4OpenGLImmediateWin32Viewer.hh
        G4OpenGLStoredWin32.hh
        G4OpenGLStoredWin32Viewer.hh
        G4OpenGLWin32Viewer.hh
    )

    list(APPEND G4VIS_MODULE_OPENGL_SOURCES
        G4OpenGLImmediateWin32.cc
        G4OpenGLImmediateWin32Viewer.cc
        G4OpenGLStoredWin32.cc
        G4OpenGLStoredWin32Viewer.cc
        G4OpenGLWin32Viewer.cc
    )

    # Add the compile definitions
    add_definitions(-DG4VIS_BUILD_OPENGLWIN32_DRIVER)

    # That should be it for Win32...
endif()



#----------------------------------------------------------------------------
# Define the Geant4 Module.
#
GEANT4_DEFINE_MODULE(NAME G4OpenGL
    HEADERS
        ${G4VIS_MODULE_OPENGL_HEADERS}
    SOURCES
        ${G4VIS_MODULE_OPENGL_SOURCES}
    GRANULAR_DEPENDENCIES
        G4UIbasic
        G4UIcommon
        G4csg
        G4event
        G4run
        G4particles
        G4processes
        G4track
        G4materials
        G4geometrymng
        G4gl2ps
        G4globman
        G4graphics_reps
        G4digits
        G4hits
        G4intercoms
        G4modeling
        G4specsolids
        G4tracking
        G4vis_management
    GLOBAL_DEPENDENCIES
        G4digits_hits
        G4event
        G4run
        G4particles
        G4processes
        G4track
        G4materials
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
        ${G4VIS_MODULE_OPENGL_LINK_LIBRARIES}
)

# List any source specific properties here

