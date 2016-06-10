#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_deex_evaporation
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_deex.G4hadronic_deex_evaporation
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 88406 2015-02-18 09:13:29Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPNumerics/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/barions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/ions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/mesons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/leptons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/fission/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/fermi_breakup/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/gem_evaporation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/photon_evaporation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_deex_evaporation
    HEADERS
        G4AlphaEvaporationChannel.hh
        G4AlphaEvaporationProbability.hh
        G4DeuteronEvaporationChannel.hh
        G4DeuteronEvaporationProbability.hh
        G4Evaporation.hh
        G4EvaporationChannel.hh
        G4EvaporationDefaultGEMFactory.hh
        G4EvaporationFactory.hh
        G4EvaporationProbability.hh
        G4He3EvaporationChannel.hh
        G4He3EvaporationProbability.hh
        G4NeutronEvaporationChannel.hh
        G4NeutronEvaporationProbability.hh
        G4ProtonEvaporationChannel.hh
        G4ProtonEvaporationProbability.hh
        G4TritonEvaporationChannel.hh
        G4TritonEvaporationProbability.hh
        G4UnstableFragmentBreakUp.hh
        G4VEvaporation.hh
    SOURCES
        G4AlphaEvaporationChannel.cc
        G4AlphaEvaporationProbability.cc
        G4DeuteronEvaporationChannel.cc
        G4DeuteronEvaporationProbability.cc
        G4Evaporation.cc
        G4EvaporationChannel.cc
        G4EvaporationDefaultGEMFactory.cc
        G4EvaporationFactory.cc
        G4EvaporationProbability.cc
        G4He3EvaporationChannel.cc
        G4He3EvaporationProbability.cc
        G4NeutronEvaporationChannel.cc
        G4NeutronEvaporationProbability.cc
        G4ProtonEvaporationChannel.cc
        G4ProtonEvaporationProbability.cc
        G4TritonEvaporationChannel.cc
        G4TritonEvaporationProbability.cc
        G4UnstableFragmentBreakUp.cc
        G4VEvaporation.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4globman
        G4had_mod_util
        G4hadronic_deex_fission
        G4hadronic_deex_gem_evaporation
        G4hadronic_deex_management
        G4hadronic_deex_photon_evaporation
        G4hadronic_deex_util
        G4hadronic_mgt
        G4hadronic_util
        G4hepnumerics
        G4ions
        G4leptons
        G4materials
        G4mesons
        G4partman
        G4track
    GLOBAL_DEPENDENCIES
        G4global
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

