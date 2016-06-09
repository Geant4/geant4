#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_deex_photon_evaporation
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_deex.G4hadronic_deex_photon_evaporation
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake,v 1.2 2010-11-29 17:48:07 bmorgan Exp $
# GEANT4 Tag $Name: not supported by cvs2svn $
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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/lowenergy/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/cross_sections/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/processes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_deex_photon_evaporation
    HEADERS
        G4ContinuumGammaDeexcitation.hh
        G4ContinuumGammaTransition.hh
        G4DiscreteGammaDeexcitation.hh
        G4DiscreteGammaTransition.hh
        G4E1Probability.hh
        G4E1SingleProbability1.hh
        G4LevelManager.hh
        G4LevelReader.hh
        G4NeutronRadCapture.hh
        G4NuclearLevel.hh
        G4NuclearLevelManager.hh
        G4NuclearLevelStore.hh
        G4NucLevel.hh
        G4PhotonEvaporation.hh
        G4PromptPhotonEvaporation.hh
        G4PtrLevelVector.hh
        G4RandGeneralTmp.hh
        G4VGammaDeexcitation.hh
        G4VGammaTransition.hh
        G4VPhotonEvaporation.hh
    SOURCES
        G4ContinuumGammaDeexcitation.cc
        G4ContinuumGammaTransition.cc
        G4DiscreteGammaDeexcitation.cc
        G4DiscreteGammaTransition.cc
        G4E1Probability.cc
        G4E1SingleProbability1.cc
        G4LevelManager.cc
        G4LevelReader.cc
        G4NeutronRadCapture.cc
        G4NuclearLevel.cc
        G4NuclearLevelManager.cc
        G4NuclearLevelStore.cc
        G4NucLevel.cc
        G4PhotonEvaporation.cc
        G4PromptPhotonEvaporation.cc
        G4VGammaDeexcitation.cc
        G4VGammaTransition.cc
        G4VPhotonEvaporation.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4emlowenergy
        G4geometrymng
        G4globman
        G4had_mod_man
        G4had_mod_util
        G4hadronic_deex_management
        G4hadronic_deex_util
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

