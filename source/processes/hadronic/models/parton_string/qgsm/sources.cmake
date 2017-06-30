#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_qgstring
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4had_string.G4hadronic_qgstring
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 102316 2017-01-20 16:12:52Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
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
include_directories(${CMAKE_SOURCE_DIR}/source/particles/shortlived/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/cross_sections/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/parton_string/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/processes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_qgstring
    HEADERS
        G4ASCCrossSection.hh
        G4AnnihilationCrossSection.hh
        G4BaryonSplitter.hh
        G4DiffractiveStringBuilder.hh
        G4GammaAnnCrossSection.hh
        G4GammaParticipants.hh
        G4MesonSplitter.hh
        G4PartonPair.hh
        G4QGSDiffractiveExcitation.hh
        G4QGSMParameters.hh
        G4QGSMSplitableHadron.hh
        G4QGSModel.hh
        G4QGSModel.icc
        G4QGSParticipants.hh
        G4QuarkExchange.hh
        G4Reggeons.hh
        G4SPBaryon.hh
        G4SPBaryonTable.hh
        G4SingleDiffractiveExcitation.hh
        G4SoftStringBuilder.hh
        G4SPPartonInfo.hh
        G4VAnnihilationCrossSection.hh
    SOURCES
        G4ASCCrossSection.cc
        G4AnnihilationCrossSection.cc
        G4BaryonSplitter.cc
        G4DiffractiveStringBuilder.cc
        G4GammaAnnCrossSection.cc
        G4GammaParticipants.cc
        G4MesonSplitter.cc
        G4PartonPair.cc
        G4QGSDiffractiveExcitation.cc
        G4QGSMParameters.cc
        G4QGSMSplitableHadron.cc
        G4QGSParticipants.cc
        G4QuarkExchange.cc
        G4Reggeons.cc
        G4SingleDiffractiveExcitation.cc
        G4SoftStringBuilder.cc
        G4SPBaryon.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4had_mod_man
        G4had_mod_util
        G4had_string_man
        G4hadronic_mgt
        G4hadronic_proc
        G4hadronic_util
        G4hadronic_xsect
        G4hepnumerics
        G4ions
        G4leptons
        G4materials
        G4mesons
        G4partman
        G4procman
        G4shortlived
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

