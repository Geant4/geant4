#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_hetcpp_utils
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4had_hetcpp.G4hadronic_hetcpp_utils
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 106962 2017-10-31 08:37:31Z gcosmo $
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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/binary_cascade/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/im_r_matrix/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/pre_equilibrium/exciton_model/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_inclxx_utils
    HEADERS
        G4INCLAllocationPool.hh
        G4INCLBook.hh
        G4INCLConfig.hh
        G4INCLConfigEnums.hh
        G4INCLEventInfo.hh
        G4INCLFinalState.hh
        G4INCLGeant4Compat.hh
        G4INCLGeant4Random.hh
        G4INCLGlobalInfo.hh
        G4INCLGlobals.hh
        G4INCLHashing.hh
        G4INCLHFB.hh
        G4INCLHornerFormEvaluator.hh
        G4INCLIAvatar.hh
        G4INCLIChannel.hh
        G4INCLIFunction1D.hh
        G4INCLInterpolationTable.hh
        G4INCLIntersection.hh
        G4INCLInvFInterpolationTable.hh
        G4INCLIRandomGenerator.hh
        G4INCLLogger.hh
        G4INCLNaturalIsotopicDistributions.hh
        G4INCLNuclearMassTable.hh
        G4INCLParticle.hh
        G4INCLParticleSpecies.hh
        G4INCLParticleTable.hh
        G4INCLParticleType.hh
        G4INCLRandom.hh
        G4INCLRandomSeedVector.hh
        G4INCLRanecu.hh
        G4INCLRanecu3.hh
        G4INCLRootFinder.hh
        G4INCLThreeVector.hh
        G4INCLUnorderedVector.hh
        G4INCLVersion.hh

    SOURCES
        G4INCLConfig.cc
        G4INCLConfigVersion.cc
        G4INCLEventInfo.cc
        G4INCLFinalState.cc
        G4INCLGlobals.cc
        G4INCLHFB.cc
        G4INCLIAvatar.cc
        G4INCLIChannel.cc
        G4INCLIFunction1D.cc
        G4INCLInterpolationTable.cc
        G4INCLInvFInterpolationTable.cc
        G4INCLLogger.cc
        G4INCLNaturalIsotopicDistributions.cc
        G4INCLNuclearMassTable.cc
        G4INCLParticle.cc
        G4INCLParticleSpecies.cc
        G4INCLParticleTable.cc
        G4INCLRandom.cc
        G4INCLRandomSeedVector.cc
        G4INCLRanecu.cc
        G4INCLRanecu3.cc
        G4INCLRootFinder.cc

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
        G4intercoms
        G4had_mod_man
        G4had_preequ_exciton
        G4hadronic_mgt
        G4hadronic_util
        G4hadronic_xsect
        G4hadronic_deex_evaporation
        G4hadronic_deex_fermi_breakup
        G4hadronic_deex_handler
        G4hadronic_deex_management
        G4hadronic_deex_multifragmentation
        G4hadronic_deex_photon_evaporation
        G4hadronic_deex_util

    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4materials
        G4particles
        G4track
        G4intercoms

LINK_LIBRARIES)# List any source specific properties here