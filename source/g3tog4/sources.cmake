#------------------------------------------------------------------------------
# sources.cmake
# Module : G3toG4
# Package: Geant4.src.G3toG4
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 88189 2015-02-02 17:22:55Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${USOLIDS_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/detector/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/magneticfield/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/Boolean/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/specific/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/decay/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)
include_directories(${CMAKE_SOURCE_DIR}/source/tracking/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G3toG4
    HEADERS
        G3DetTable.hh
        G3DetTableEntry.hh
        G3Division.hh
        G3EleTable.hh
        G3G4Interface.hh
        G3MatTable.hh
        G3MatTableEntry.hh
        G3MedTable.hh
        G3MedTableEntry.hh
        G3PartTable.hh
        G3Pos.hh
        G3RotTable.hh
        G3RotTableEntry.hh
        G3SensVolVector.hh
        G3VolTable.hh
        G3VolTableEntry.hh
        G3toG4.hh
        G3toG4BuildTree.hh
        G3toG4Defs.hh
        G3toG4MANY.hh
        G3toG4MakeSolid.hh
        G3toG4RotationMatrix.hh
    SOURCES
        G3DetTable.cc
        G3DetTableEntry.cc
        G3Division.cc
        G3EleTable.cc
        G3MatTable.cc
        G3MatTableEntry.cc
        G3MedTable.cc
        G3MedTableEntry.cc
        G3NegVolPars.cc
        G3PartTable.cc
        G3Pos.cc
        G3RotTable.cc
        G3RotTableEntry.cc
        G3VolTable.cc
        G3VolTableEntry.cc
        G3toG4BuildTree.cc
        G3toG4MANY.cc
        G3toG4MakeSolid.cc
        G3toG4RotationMatrix.cc
        G4BuildGeom.cc
        G4ggclos.cc
        G4gsatt.cc
        G4gsbool.cc
        G4gsdet.cc
        G4gsdeta.cc
        G4gsdetd.cc
        G4gsdeth.cc
        G4gsdetu.cc
        G4gsdetv.cc
        G4gsdk.cc
        G4gsdvn.cc
        G4gsdvn2.cc
        G4gsdvt.cc
        G4gsdvt2.cc
        G4gsdvx.cc
        G4gsmate.cc
        G4gsmixt.cc
        G4gspart.cc
        G4gspos.cc
        G4gsposp.cc
        G4gsrotm.cc
        G4gstmed.cc
        G4gstpar.cc
        G4gsvolu.cc
        clparse.cc
    GRANULAR_DEPENDENCIES
        G4brep
        G4csg
        G4decay
        G4detector
        G4geomBoolean
        G4geometrymng
        G4globman
        G4graphics_reps
        G4hits
        G4magneticfield
        G4materials
        G4partman
        G4procman
        G4specsolids
        G4track
        G4tracking
        G4volumes
    GLOBAL_DEPENDENCIES
        G4digits_hits
        G4geometry
        G4global
        G4graphics_reps
        G4materials
        G4particles
        G4processes
        G4track
        G4tracking
    LINK_LIBRARIES
)

# List any source specific properties here

