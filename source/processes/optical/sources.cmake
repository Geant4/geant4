#------------------------------------------------------------------------------
# sources.cmake
# Module : G4optical
# Package: Geant4.src.G4processes.G4optical
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 73616 2013-09-02 10:18:24Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/magneticfield/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/scoring/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/detector/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4optical
    HEADERS
        G4OpAbsorption.hh
        G4OpBoundaryProcess.hh
        G4OpMieHG.hh
        G4OpProcessSubType.hh
        G4OpRayleigh.hh
        G4OpWLS.hh
        G4VWLSTimeGeneratorProfile.hh
        G4WLSTimeGeneratorProfileDelta.hh
        G4WLSTimeGeneratorProfileExponential.hh
    SOURCES
        G4OpAbsorption.cc
        G4OpBoundaryProcess.cc
        G4OpMieHG.cc
        G4OpRayleigh.cc
        G4OpWLS.cc
        G4VWLSTimeGeneratorProfile.cc
        G4WLSTimeGeneratorProfileDelta.cc
        G4WLSTimeGeneratorProfileExponential.cc
    GRANULAR_DEPENDENCIES
        G4bosons
        G4geometrymng
        G4globman
        G4materials
        G4navigation
        G4partman
        G4procman
        G4track
        G4volumes
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

