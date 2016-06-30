#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_deex_util
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_deex.G4hadronic_deex_util
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 96634 2016-04-27 09:31:49Z gcosmo $
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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_deex_util
    HEADERS
        G4AlphaCoulombBarrier.hh
        G4CameronGilbertPairingCorrections.hh
        G4CameronGilbertShellCorrections.hh
        G4CameronShellPlusPairingCorrections.hh
        G4CameronTruranHilfPairingCorrections.hh
        G4CameronTruranHilfShellCorrections.hh
	G4ChatterjeeCrossSection.hh
        G4ConstantLevelDensityParameter.hh
        G4CookPairingCorrections.hh
        G4CookShellCorrections.hh
        G4CoulombBarrier.hh
        G4DeuteronCoulombBarrier.hh
        G4EvaporationLevelDensityParameter.hh
        G4He3CoulombBarrier.hh
        G4KalbachCrossSection.hh
        G4NeutronCoulombBarrier.hh
        G4PairingCorrection.hh
        G4ProtonCoulombBarrier.hh
        G4ShellCorrection.hh
        G4TritonCoulombBarrier.hh
        G4VCoulombBarrier.hh
        G4VLevelDensityParameter.hh
    SOURCES
        G4AlphaCoulombBarrier.cc
        G4CameronGilbertPairingCorrections.cc
        G4CameronGilbertShellCorrections.cc
        G4CameronShellPlusPairingCorrections.cc
        G4CameronTruranHilfPairingCorrections.cc
        G4CameronTruranHilfShellCorrections.cc
        G4ConstantLevelDensityParameter.cc
        G4CookPairingCorrections.cc
        G4CookShellCorrections.cc
        G4CoulombBarrier.cc
        G4DeuteronCoulombBarrier.cc
        G4EvaporationLevelDensityParameter.cc
        G4He3CoulombBarrier.cc
        G4NeutronCoulombBarrier.cc
        G4PairingCorrection.cc
        G4ProtonCoulombBarrier.cc
        G4ShellCorrection.cc
        G4TritonCoulombBarrier.cc
        G4VCoulombBarrier.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4globman
        G4had_mod_util
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

