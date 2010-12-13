#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_HE
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_HE
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake,v 1.2 2010-12-13 12:08:19 bmorgan Exp $
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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_HE
    HEADERS
        G4HEAntiKaonZeroInelastic.hh
        G4HEAntiLambdaInelastic.hh
        G4HEAntiNeutronInelastic.hh
        G4HEAntiOmegaMinusInelastic.hh
        G4HEAntiProtonInelastic.hh
        G4HEAntiSigmaMinusInelastic.hh
        G4HEAntiSigmaPlusInelastic.hh
        G4HEAntiSigmaZeroInelastic.hh
        G4HEAntiXiMinusInelastic.hh
        G4HEAntiXiZeroInelastic.hh
        G4HEInelastic.hh
        G4HEKaonMinusInelastic.hh
        G4HEKaonPlusInelastic.hh
        G4HEKaonZeroInelastic.hh
        G4HEKaonZeroLongInelastic.hh
        G4HEKaonZeroShortInelastic.hh
        G4HELambdaInelastic.hh
        G4HENeutronInelastic.hh
        G4HEOmegaMinusInelastic.hh
        G4HEPionMinusInelastic.hh
        G4HEPionPlusInelastic.hh
        G4HEProtonInelastic.hh
        G4HESigmaMinusInelastic.hh
        G4HESigmaPlusInelastic.hh
        G4HESigmaZeroInelastic.hh
        G4HEVector.hh
        G4HEXiMinusInelastic.hh
        G4HEXiZeroInelastic.hh
    SOURCES
        G4HEAntiKaonZeroInelastic.cc
        G4HEAntiLambdaInelastic.cc
        G4HEAntiNeutronInelastic.cc
        G4HEAntiOmegaMinusInelastic.cc
        G4HEAntiProtonInelastic.cc
        G4HEAntiSigmaMinusInelastic.cc
        G4HEAntiSigmaPlusInelastic.cc
        G4HEAntiSigmaZeroInelastic.cc
        G4HEAntiXiMinusInelastic.cc
        G4HEAntiXiZeroInelastic.cc
        G4HEInelastic.cc
        G4HEKaonMinusInelastic.cc
        G4HEKaonPlusInelastic.cc
        G4HEKaonZeroInelastic.cc
        G4HEKaonZeroLongInelastic.cc
        G4HEKaonZeroShortInelastic.cc
        G4HELambdaInelastic.cc
        G4HENeutronInelastic.cc
        G4HEOmegaMinusInelastic.cc
        G4HEPionMinusInelastic.cc
        G4HEPionPlusInelastic.cc
        G4HEProtonInelastic.cc
        G4HESigmaMinusInelastic.cc
        G4HESigmaPlusInelastic.cc
        G4HESigmaZeroInelastic.cc
        G4HEVector.cc
        G4HEXiMinusInelastic.cc
        G4HEXiZeroInelastic.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4had_mod_man
        G4hadronic_mgt
        G4hadronic_util
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

