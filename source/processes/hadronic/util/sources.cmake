#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_util
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_util
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 103898 2017-05-03 08:25:03Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/barions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/ions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/mesons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/leptons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_util
    HEADERS
        DumpFrame.hh
        G4Bessel.hh
        G4Delete.hh
        G4DynamicParticleVector.hh
        G4GHEKinematicsVector.hh
        G4HadFinalState.hh
        G4HadParticleCodes.hh
        G4HadProjectile.hh
        G4HadReentrentException.hh
        G4HadSecondary.hh
        G4HadSignalHandler.hh
        G4HadTmpUtil.hh
	G4HadronicDeprecate.hh
	G4HadronicDeveloperParameters.hh
        G4HadronicException.hh
        G4HadronicWhiteBoard.hh
        G4IsoResult.hh
        G4LightMedia.hh
        G4NuclearPolarization.hh
        G4Nucleus.hh
        G4Pair.hh
        G4ReactionProduct.hh
        G4ReactionProductVector.hh
        G4StableIsotopes.hh
        G4ping.hh
    SOURCES
        G4Bessel.cc
        G4HadFinalState.cc
        G4HadProjectile.cc
        G4HadSecondary.cc
        G4HadSignalHandler.cc
        G4HadTmpUtil.cc
	G4HadronicDeveloperParameters.cc
        G4HadronicWhiteBoard.cc
        G4IsoResult.cc
        G4LightMedia.cc
        G4NuclearPolarization.cc
        G4Nucleus.cc
        G4ReactionProduct.cc
        G4StableIsotopes.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4ions
        G4leptons
        G4materials
        G4mesons
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

