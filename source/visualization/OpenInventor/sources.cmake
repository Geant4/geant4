#------------------------------------------------------------------------------
# sources.cmake
# Module : G4OpenInventor
# Package: Geant4.src.G4visualization.G4OpenInventor
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 102233 2017-01-13 15:59:06Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${USOLIDS_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/specific/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/interfaces/common/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/tracking/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/externals/gl2ps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/modeling/include)


#----------------------------------------------------------------------------
# Generic Inventor Headers, base library, OpenGL and Geant4 defines
#
include_directories(${INVENTOR_INCLUDE_DIR})
set(G4VIS_MODULE_OPENINVENTOR_LINK_LIBRARIES ${INVENTOR_LIBRARY})

include_directories(${OPENGL_INCLUDE_DIR})
list(APPEND G4VIS_MODULE_OPENINVENTOR_LINK_LIBRARIES ${OPENGL_LIBRARIES})

add_definitions(-DG4VIS_BUILD_OI_DRIVER)



#----------------------------------------------------------------------------
# Geant4 OpenInventor Core sources and headers (all platforms)
#
set(G4VIS_MODULE_OPENINVENTOR_HEADERS
  G4OpenInventor.hh
  G4OpenInventorSceneHandler.hh
  G4OpenInventorTransform3D.hh
  G4OpenInventorViewer.hh
  G4VisFeaturesOfOpenInventor.hh
  Geant4_SoPolyhedron.h
  SoG4LineSet.h
  SoG4MarkerSet.h
  SoG4Polyhedron.h
  )

set(G4VIS_MODULE_OPENINVENTOR_SOURCES
  G4OpenInventor.cc
  G4OpenInventorSceneHandler.cc
  G4OpenInventorTransform3D.cc
  G4OpenInventorViewer.cc
  G4VisFeaturesOfOpenInventor.cc
  SbPainter.cc
  SbPainterPS.cc
  SoAlternateRepAction.cc
  SoBox.cc
  SoCons.cc
  SoCounterAction.cc
  SoDetectorTreeKit.cc
  SoGL2PSAction.cc
  SoImageWriter.cc
  SoMarkerSet.cc
  SoPolyhedron.cc
  SoStyleCache.cc
  SoTrap.cc
  SoTrd.cc
  SoTubs.cc
  )



#----------------------------------------------------------------------------
# UNIX Only (Xt) sources
#
if(UNIX)
  list(APPEND G4VIS_MODULE_OPENINVENTOR_HEADERS
    G4OpenInventorX.hh
    G4OpenInventorXt.hh
    G4OpenInventorXtExaminerViewerMessenger.hh
    G4OpenInventorXtExaminerViewer.hh
    G4OpenInventorXtExtended.hh
    G4OpenInventorXtExtendedViewer.hh
    G4OpenInventorXtViewer.hh
    wheelmouse.h
    SoXtInternal.h
    console.h
    favorites.h
    saveViewPt.h
    pickext.h
    pickref.h
    wireframe.h
    )

  list(APPEND G4VIS_MODULE_OPENINVENTOR_SOURCES
    G4OpenInventorXt.cc
    G4OpenInventorXtExaminerViewer.cc
    G4OpenInventorXtExaminerViewerMessenger.cc
    G4OpenInventorXtExtended.cc
    G4OpenInventorXtExtendedViewer.cc
    G4OpenInventorXtViewer.cc
    wheelmouse.cc
    )

  # Add the definitions for SoXt
  add_definitions(-DG4INTY_BUILD_XT)
  add_definitions(-DG4VIS_BUILD_OIX_DRIVER)

  # SoXt Library
  list(APPEND G4VIS_MODULE_OPENINVENTOR_LINK_LIBRARIES
    ${INVENTOR_SOXT_LIBRARY}
    )

  # We also need Xm and X11
  include_directories(${X11_INCLUDE_DIR})
  include_directories(${MOTIF_INCLUDE_DIR})
  list(APPEND G4VIS_MODULE_OPENINVENTOR_LINK_LIBRARIES
   ${MOTIF_LIBRARIES}
   ${X11_LIBRARIES}
   ${X11_Xpm_LIB}
   )
endif()


#----------------------------------------------------------------------------
# WIN32 Only (Win32) sources
#
if(WIN32)
  set(G4VIS_MODULE_OPENINVENTOR_HEADERS
    ${G4VIS_MODULE_OPENINVENTOR_HEADERS}
    G4OpenInventorWin.hh
    G4OpenInventorWin32.hh
    G4OpenInventorWinViewer.hh
    )

  set(G4VIS_MODULE_OPENINVENTOR_SOURCES
    ${G4VIS_MODULE_OPENINVENTOR_SOURCES}
    G4OpenInventorWin.cc
    G4OpenInventorWinViewer.cc
    )

  # Add the include for SoWin

  # Add the definitions for SoWin
  add_definitions(-DG4INTY_BUILD_WIN32)
  add_definitions(-DG4VIS_BUILD_OIWIN32_DRIVER)

  # SoWin Library
  list(APPEND G4VIS_MODULE_OPENINVENTOR_LINK_LIBRARIES
    ${INVENTOR_SOWIN_LIBRARY}
    )
endif()



#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4OpenInventor
  HEADERS
    ${G4VIS_MODULE_OPENINVENTOR_HEADERS}
  SOURCES
    ${G4VIS_MODULE_OPENINVENTOR_SOURCES}
  GRANULAR_DEPENDENCIES
    G4UIcommon
    G4csg
    G4geometrymng
    G4gl2ps
    G4globman
    G4graphics_reps
    G4hits
    G4intercoms
    G4materials
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
    G4materials
    G4modeling
    G4tracking
    G4vis_management
  LINK_LIBRARIES
    ${G4VIS_MODULE_OPENINVENTOR_LINK_LIBRARIES}
  )

# List any source specific properties here

