#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_LE
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_LE
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake,v 1.1 2010-09-29 19:03:39 bmorgan Exp $
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
GEANT4_DEFINE_MODULE(NAME G4hadronic_LE
    HEADERS
        G4LCapture.hh
        G4LEAlphaInelastic.hh
        G4LEAntiKaonZeroInelastic.hh
        G4LEAntiLambdaInelastic.hh
        G4LEAntiNeutronInelastic.hh
        G4LEAntiOmegaMinusInelastic.hh
        G4LEAntiProtonInelastic.hh
        G4LEAntiSigmaMinusInelastic.hh
        G4LEAntiSigmaPlusInelastic.hh
        G4LEAntiXiMinusInelastic.hh
        G4LEAntiXiZeroInelastic.hh
        G4LEDeuteronInelastic.hh
        G4LEKaonMinusInelastic.hh
        G4LEKaonPlusInelastic.hh
        G4LEKaonZeroInelastic.hh
        G4LEKaonZeroLInelastic.hh
        G4LEKaonZeroSInelastic.hh
        G4LELambdaInelastic.hh
        G4LENeutronInelastic.hh
        G4LEOmegaMinusInelastic.hh
        G4LEPionMinusInelastic.hh
        G4LEPionPlusInelastic.hh
        G4LEProtonInelastic.hh
        G4LESigmaMinusInelastic.hh
        G4LESigmaPlusInelastic.hh
        G4LETritonInelastic.hh
        G4LEXiMinusInelastic.hh
        G4LEXiZeroInelastic.hh
        G4LElastic.hh
        G4LFission.hh
    SOURCES
        G4LCapture.cc
        G4LEAlphaInelastic.cc
        G4LEAntiKaonZeroInelastic.cc
        G4LEAntiLambdaInelastic.cc
        G4LEAntiNeutronInelastic.cc
        G4LEAntiOmegaMinusInelastic.cc
        G4LEAntiProtonInelastic.cc
        G4LEAntiSigmaMinusInelastic.cc
        G4LEAntiSigmaPlusInelastic.cc
        G4LEAntiXiMinusInelastic.cc
        G4LEAntiXiZeroInelastic.cc
        G4LEDeuteronInelastic.cc
        G4LEKaonMinusInelastic.cc
        G4LEKaonPlusInelastic.cc
        G4LEKaonZeroInelastic.cc
        G4LELambdaInelastic.cc
        G4LENeutronInelastic.cc
        G4LEOmegaMinusInelastic.cc
        G4LEPionMinusInelastic.cc
        G4LEPionPlusInelastic.cc
        G4LEProtonInelastic.cc
        G4LESigmaMinusInelastic.cc
        G4LESigmaPlusInelastic.cc
        G4LETritonInelastic.cc
        G4LEXiMinusInelastic.cc
        G4LEXiZeroInelastic.cc
        G4LElastic.cc
        G4LFission.cc
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

