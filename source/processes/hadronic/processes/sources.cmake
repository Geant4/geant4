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
# $Id: sources.cmake,v 1.1 2010-09-29 19:08:07 bmorgan Exp $
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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/cross_sections/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_proc
    HEADERS
        G4AlphaInelasticProcess.hh
        G4AntiLambdaInelasticProcess.hh
        G4AntiNeutronInelasticProcess.hh
        G4AntiOmegaMinusInelasticProcess.hh
        G4AntiProtonInelasticProcess.hh
        G4AntiSigmaMinusInelasticProcess.hh
        G4AntiSigmaPlusInelasticProcess.hh
        G4AntiXiMinusInelasticProcess.hh
        G4AntiXiZeroInelasticProcess.hh
        G4DeuteronInelasticProcess.hh
        G4ElectronNuclearProcess.hh
        G4HadronCaptureProcess.hh
        G4HadronElasticProcess.hh
        G4HadronFissionProcess.hh
        G4IonInelasticProcess.hh
        G4KaonMinusInelasticProcess.hh
        G4KaonPlusInelasticProcess.hh
        G4KaonZeroLInelasticProcess.hh
        G4KaonZeroSInelasticProcess.hh
        G4LambdaInelasticProcess.hh
        G4NeutronInelasticProcess.hh
        G4OmegaMinusInelasticProcess.hh
        G4PhotoNuclearProcess.hh
        G4PionMinusInelasticProcess.hh
        G4PionPlusInelasticProcess.hh
        G4PositronNuclearProcess.hh
        G4ProtonInelasticProcess.hh
        G4SigmaMinusInelasticProcess.hh
        G4SigmaPlusInelasticProcess.hh
        G4TritonInelasticProcess.hh
        G4XiMinusInelasticProcess.hh
        G4XiZeroInelasticProcess.hh
    SOURCES
        G4ElectronNuclearProcess.cc
        G4HadronCaptureProcess.cc
        G4HadronElasticProcess.cc
        G4HadronFissionProcess.cc
        G4PhotoNuclearProcess.cc
        G4PositronNuclearProcess.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4had_mod_man
        G4hadronic_mgt
        G4hadronic_util
        G4hadronic_xsect
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

