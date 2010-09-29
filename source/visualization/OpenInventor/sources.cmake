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
# $Id: sources.cmake,v 1.1 2010-09-29 19:12:30 bmorgan Exp $
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

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4OpenInventor
    HEADERS
        G4OpenInventor.hh
        G4OpenInventorSceneHandler.hh
        G4OpenInventorTransform3D.hh
        G4OpenInventorViewer.hh
        G4OpenInventorWin.hh
        G4OpenInventorWin32.hh
        G4OpenInventorWinViewer.hh
        G4OpenInventorX.hh
        G4OpenInventorXt.hh
        G4OpenInventorXtViewer.hh
        Geant4_SoPolyhedron.h
        SoG4LineSet.h
        SoG4MarkerSet.h
        SoG4Polyhedron.h
    SOURCES
        G4OpenInventor.cc
        G4OpenInventorSceneHandler.cc
        G4OpenInventorTransform3D.cc
        G4OpenInventorViewer.cc
        G4OpenInventorWin.cc
        G4OpenInventorWinViewer.cc
        G4OpenInventorXt.cc
        G4OpenInventorXtViewer.cc
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
)

# List any source specific properties here

