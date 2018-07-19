#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_deex_fermi_breakup
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_deex.G4hadronic_deex_fermi_breakup
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 102988 2017-03-07 16:22:20Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPNumerics/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/barions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/leptons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_deex_fermi_breakup
    HEADERS
        G4FermiBreakUpVI.hh
        G4FermiChannels.hh
        G4FermiDecayProbability.hh
        G4FermiFragment.hh
        G4FermiFragmentsPoolVI.hh
        G4FermiPair.hh
        G4FermiPhaseSpaceDecay.hh
        G4VFermiBreakUp.hh
    SOURCES
        G4FermiBreakUpVI.cc
        G4FermiDecayProbability.cc
        G4FermiFragment.cc
        G4FermiFragmentsPoolVI.cc
        G4FermiPair.cc
        G4FermiPhaseSpaceDecay.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4globman
        G4had_mod_util
        G4hadronic_deex_management
        G4hadronic_util
        G4hepnumerics
        G4leptons
        G4partman
    GLOBAL_DEPENDENCIES
        G4global
        G4particles
    LINK_LIBRARIES
)

# List any source specific properties here

