#------------------------------------------------------------------------------
# sources.cmake
# Module : G4geomdivision
# Package: Geant4.src.G4geometry.G4geomdivision
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 79490 2014-03-05 15:23:33Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/specific/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4geomdivision
    HEADERS
        G4PVDivision.hh
        G4PVDivisionFactory.hh
        G4ParameterisationBox.hh
        G4ParameterisationCons.hh
        G4ParameterisationPara.hh
        G4ParameterisationPolycone.hh
        G4ParameterisationPolyhedra.hh
        G4ParameterisationTrd.hh
        G4ParameterisationTubs.hh
        G4ReplicatedSlice.hh
        G4VDivisionParameterisation.hh
        G4VDivisionParameterisation.icc
    SOURCES
        G4PVDivision.cc
        G4PVDivisionFactory.cc
        G4ParameterisationBox.cc
        G4ParameterisationCons.cc
        G4ParameterisationPara.cc
        G4ParameterisationPolycone.cc
        G4ParameterisationPolyhedra.cc
        G4ParameterisationTrd.cc
        G4ParameterisationTubs.cc
        G4ReplicatedSlice.cc
        G4VDivisionParameterisation.cc
    GRANULAR_DEPENDENCIES
        G4csg
        G4geometrymng
        G4globman
        G4specsolids
    GLOBAL_DEPENDENCIES
        G4global
    LINK_LIBRARIES
)

# List any source specific properties here

