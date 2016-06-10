#-----------------------------------------------------------------------
# sources.cmake
# Module : G4visHepRep
# Package: Geant4.src.G4visualization.G4visHepRep
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 88190 2015-02-02 17:24:54Z gcosmo $
#
#-----------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${ZLIB_INCLUDE_DIRS})
include_directories(${USOLIDS_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/specific/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/tracking/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/modeling/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4visHepRep
  HEADERS
    G4HepRep.hh
    G4HepRepFile.hh
    G4HepRepFileSceneHandler.hh
    G4HepRepFileViewer.hh
    G4HepRepFileXMLWriter.hh
    G4HepRepMessenger.hh
    G4HepRepSceneHandler.hh
    G4HepRepViewer.hh
  SOURCES
    BHepRepWriter.cc
    DefaultHepRep.cc
    DefaultHepRepAction.cc
    DefaultHepRepAttDef.cc
    DefaultHepRepAttValue.cc
    DefaultHepRepAttribute.cc
    DefaultHepRepDefinition.cc
    DefaultHepRepFactory.cc
    DefaultHepRepInstance.cc
    DefaultHepRepInstanceTree.cc
    DefaultHepRepPoint.cc
    DefaultHepRepTreeID.cc
    DefaultHepRepType.cc
    DefaultHepRepTypeTree.cc
    DeflateOutputStreamBuffer.cc
    G4HepRep.cc
    G4HepRepFile.cc
    G4HepRepFileSceneHandler.cc
    G4HepRepFileViewer.cc
    G4HepRepFileXMLWriter.cc
    G4HepRepMessenger.cc
    G4HepRepSceneHandler.cc
    G4HepRepViewer.cc
    GZIPOutputStream.cc
    GZIPOutputStreamBuffer.cc
    IndentPrintWriter.cc
    XMLHepRepFactory.cc
    XMLHepRepWriter.cc
    XMLWriter.cc
    ZipOutputStream.cc
    ZipOutputStreamBuffer.cc
  GRANULAR_DEPENDENCIES
    G4csg
    G4geometrymng
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
    G4global
    G4graphics_reps
    G4intercoms
    G4materials
    G4modeling
    G4tracking
    G4vis_management
  LINK_LIBRARIES
    ${ZLIB_LIBRARIES}
  )

# List any source specific properties here

