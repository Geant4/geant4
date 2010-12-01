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
# $Id: sources.cmake,v 1.5 2010-12-01 16:57:53 bmorgan Exp $
# GEANT4 Tag $Name: not supported by cvs2svn $
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
    G4OpenGLBitMapStore.hh
    G4OpenGLFontBaseStore.hh
    G4OpenGLSceneHandler.hh
    G4OpenGLSceneHandler.icc
    G4OpenGLViewerMessenger.hh
    G4OpenGLTransform3D.hh
)

set(G4VIS_MODULE_OPENGL_SOURCES 
    G4OpenGLImmediateViewer.cc
    G4OpenGLImmediateSceneHandler.cc
    G4OpenGLViewer.cc
    G4OpenGLStoredViewer.cc
    G4OpenGLStoredSceneHandler.cc
    G4OpenGLBitMapStore.cc
    G4OpenGLFontBaseStore.cc
    G4OpenGLSceneHandler.cc
    G4OpenGLViewerMessenger.cc
    G4OpenGLTransform3D.cc
)

#
# May need OpenGL include here
#
include_directories(${OPENGL_INCLUDE_DIR})
set(G4VIS_MODULE_OPENGL_LINK_LIBRARIES ${OPENGL_LIBRARIES})

#
# All of the above files must have the G4VIS_BUILD_OPENGL_DRIVER definition
#
GEANT4_ADD_COMPILE_DEFINITIONS(SOURCES ${G4VIS_MODULE_OPENGL_SOURCES}
    COMPILE_DEFINITIONS G4VIS_BUILD_OPENGL_DRIVER)


#
# Qt only if selected
#
if(GEANT4_USE_QT)
    #
    # Add in the extra Qt GL sources
    #
    list(APPEND G4VIS_MODULE_OPENGL_HEADERS
        G4OpenGLImmediateQt.hh
        G4OpenGLImmediateQtViewer.hh
        G4OpenGLQtExportDialog.hh
        G4OpenGLQtMovieDialog.hh
        G4OpenGLQtViewer.hh
        G4OpenGLStoredQt.hh
        G4OpenGLStoredQtViewer.hh)

    list(APPEND G4VIS_MODULE_OPENGL_SOURCES
        G4OpenGLImmediateQt.cc
        G4OpenGLImmediateQtViewer.cc
        G4OpenGLQtExportDialog.cc
        G4OpenGLQtMovieDialog.cc
        G4OpenGLQtViewer.cc
        G4OpenGLStoredQt.cc
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

    # Add the definitions
    # Argh.. Have to remember about INTY and UI because of their use...
    GEANT4_ADD_COMPILE_DEFINITIONS(SOURCES ${G4VIS_MODULE_OPENGL_SOURCES}
        COMPILE_DEFINITIONS
        G4VIS_BUILD_OPENGL_DRIVER;G4VIS_BUILD_OPENGLQT_DRIVER;G4INTY_BUILD_QT;G4UI_BUILD_QT_SESSION)

    # And for the moc files because these are in the build tree
    set_source_files_properties(${G4OPENGL_MOC_SOURCES}
       PROPERTIES 
       COMPILE_DEFINITIONS
       "G4VIS_BUILD_OPENGL_DRIVER;G4VIS_BUILD_OPENGLQT_DRIVER"
    )

    # Add in Qt libraries
    list(APPEND G4VIS_MODULE_OPENGL_LINK_LIBRARIES ${QT_LIBRARIES})

endif()

#
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
        ${G4VIS_MODULE_OPENGL_LINK_LIBRARIES}
)

# List any source specific properties here

