#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_proc
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_proc
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 93608 2015-10-27 08:50:25Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/magneticfield/include)
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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/cuts/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/cross_sections/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/scoring/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/detector/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_proc
    HEADERS
        G4AlphaInelasticProcess.hh
        G4AntiAlphaInelasticProcess.hh
        G4AntiDeuteronInelasticProcess.hh
        G4AntiHe3InelasticProcess.hh
        G4AntiLambdaInelasticProcess.hh
        G4AntiNeutronInelasticProcess.hh
        G4AntiOmegaMinusInelasticProcess.hh
        G4AntiProtonInelasticProcess.hh
        G4AntiSigmaMinusInelasticProcess.hh
        G4AntiSigmaPlusInelasticProcess.hh
        G4AntiTritonInelasticProcess.hh
        G4AntiXiMinusInelasticProcess.hh
        G4AntiXiZeroInelasticProcess.hh
        G4DeuteronInelasticProcess.hh
        G4ElectronNuclearProcess.hh
        G4HadronCaptureProcess.hh
        G4HadronElasticProcess.hh
        G4HadronFissionProcess.hh
        G4He3InelasticProcess.hh
        G4IonInelasticProcess.hh
        G4KaonMinusInelasticProcess.hh
        G4KaonPlusInelasticProcess.hh
        G4KaonZeroLInelasticProcess.hh
        G4KaonZeroSInelasticProcess.hh
        G4LambdaInelasticProcess.hh
        G4MuonNuclearProcess.hh
        G4NeutronInelasticProcess.hh
        G4OmegaMinusInelasticProcess.hh
	G4PhotoCaptureProcess.hh
	G4PhotoFissionProcess.hh
        G4PhotoNuclearProcess.hh
        G4PionMinusInelasticProcess.hh
        G4PionPlusInelasticProcess.hh
        G4PositronNuclearProcess.hh
        G4ProtonInelasticProcess.hh
        G4SigmaMinusInelasticProcess.hh
        G4SigmaPlusInelasticProcess.hh
        G4TritonInelasticProcess.hh
        G4UCNProcessSubType.hh
        G4UCNBoundaryProcess.hh
        G4UCNBoundaryProcessMessenger.hh
        G4UCNLoss.hh
        G4UCNAbsorption.hh
        G4UCNMultiScattering.hh
        G4XiMinusInelasticProcess.hh
        G4XiZeroInelasticProcess.hh
    SOURCES
        G4AlphaInelasticProcess.cc
        G4AntiAlphaInelasticProcess.cc
        G4AntiDeuteronInelasticProcess.cc
        G4AntiHe3InelasticProcess.cc
        G4AntiLambdaInelasticProcess.cc
        G4AntiNeutronInelasticProcess.cc
        G4AntiOmegaMinusInelasticProcess.cc
        G4AntiProtonInelasticProcess.cc
        G4AntiSigmaMinusInelasticProcess.cc
        G4AntiSigmaPlusInelasticProcess.cc
        G4AntiTritonInelasticProcess.cc
        G4AntiXiMinusInelasticProcess.cc
        G4AntiXiZeroInelasticProcess.cc
        G4DeuteronInelasticProcess.cc
        G4ElectronNuclearProcess.cc
        G4HadronCaptureProcess.cc
        G4HadronElasticProcess.cc
        G4HadronFissionProcess.cc
        G4He3InelasticProcess.cc
        G4IonInelasticProcess.cc
        G4KaonMinusInelasticProcess.cc
        G4KaonPlusInelasticProcess.cc
        G4KaonZeroLInelasticProcess.cc
        G4KaonZeroSInelasticProcess.cc
        G4LambdaInelasticProcess.cc
        G4MuonNuclearProcess.cc
        G4NeutronInelasticProcess.cc
        G4OmegaMinusInelasticProcess.cc
	G4PhotoCaptureProcess.cc
	G4PhotoFissionProcess.cc
        G4PhotoNuclearProcess.cc
        G4PionMinusInelasticProcess.cc
        G4PionPlusInelasticProcess.cc
        G4PositronNuclearProcess.cc
        G4ProtonInelasticProcess.cc
        G4SigmaMinusInelasticProcess.cc
        G4SigmaPlusInelasticProcess.cc
        G4TritonInelasticProcess.cc
        G4UCNBoundaryProcess.cc
        G4UCNBoundaryProcessMessenger.cc
        G4UCNLoss.cc
        G4UCNAbsorption.cc
        G4UCNMultiScattering.cc
        G4XiMinusInelasticProcess.cc
        G4XiZeroInelasticProcess.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4cuts
        G4geometrymng
        G4globman
        G4had_mod_man
        G4hadronic_mgt
        G4hadronic_util
        G4hadronic_xsect
        G4magneticfield
        G4ions
        G4leptons
        G4mesons
        G4materials
        G4navigation
        G4partman
        G4procman
        G4track
        G4volumes
        G4digits
        G4hits
        G4scoring
        G4intercoms
    GLOBAL_DEPENDENCIES
        G4digits_hits
        G4geometry
        G4global
        G4materials
        G4particles
        G4track
        G4intercoms
    LINK_LIBRARIES
)

# List any source specific properties here

