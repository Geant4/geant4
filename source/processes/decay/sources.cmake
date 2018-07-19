#------------------------------------------------------------------------------
# sources.cmake
# Module : G4decay
# Package: Geant4.src.G4processes.G4decay
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 106151 2017-09-14 06:43:04Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/magneticfield/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadrons/ions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4decay
    HEADERS
        G4Decay.hh
        G4DecayProcessType.hh
        G4DecayWithSpin.hh
        G4PionDecayMakeSpin.hh
        G4UnknownDecay.hh
        G4VExtDecayer.hh
    SOURCES
        G4Decay.cc
        G4DecayWithSpin.cc
        G4PionDecayMakeSpin.cc
        G4UnknownDecay.cc
    GRANULAR_DEPENDENCIES
        G4csg
        G4geometrymng
        G4globman
        G4intercoms
        G4magneticfield
        G4materials
        G4navigation
        G4partman
        G4procman
        G4track
        G4volumes
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4intercoms
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

