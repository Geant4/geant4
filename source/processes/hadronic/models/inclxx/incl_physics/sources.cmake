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
# $Id: sources.cmake 85988 2014-11-06 15:38:49Z gcosmo $
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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_inclxx_physics
    HEADERS
        G4INCLPauliBlocking.hh
        G4INCLINuclearPotential.hh
        G4INCLInteractionAvatar.hh
        G4INCLTransmissionChannel.hh
        G4INCLElasticChannel.hh
        G4INCLDeltaDecayChannel.hh
        G4INCLNucleus.hh
        G4INCLCrossSectionsINCL46.hh
        G4INCLCluster.hh
        G4INCLParticleEntryChannel.hh
        G4INCLNDFGaussian.hh
        G4INCLKinematicsUtils.hh
        G4INCLClustering.hh
        G4INCLNuclearDensityFactory.hh
        G4INCLSurfaceAvatar.hh
        G4INCLCascadeAction.hh
        G4INCLPhaseSpaceRauboldLynch.hh
        G4INCLReflectionChannel.hh
        G4INCLBinaryCollisionAvatar.hh
        G4INCLNDFHardSphere.hh
        G4INCLStore.hh
        G4INCLNDFWoodsSaxon.hh
        G4INCLPauliStandard.hh
        G4INCLDeuteronDensity.hh
        G4INCLPauliGlobal.hh
        G4INCLPiNElasticChannel.hh
        G4INCLClusterDecay.hh
        G4INCLProjectileRemnant.hh
        G4INCLIPauli.hh
        G4INCLIPhaseSpaceGenerator.hh
        G4INCLDecayAvatar.hh
        G4INCLPiNToMultiPionsChannel.hh
        G4INCLNuclearPotentialEnergyIsospin.hh
        G4INCLAvatarDumpAction.hh
        G4INCLCoulombDistortion.hh
        G4INCLClusteringModelNone.hh
        G4INCLPhaseSpaceKopylov.hh
        G4INCLIClusteringModel.hh
        G4INCLDeJongSpin.hh
        G4INCLParticleEntryAvatar.hh
        G4INCLCoulombNone.hh
        G4INCLIPropagationModel.hh
        G4INCLNuclearPotentialEnergyIsospinSmooth.hh
        G4INCLStandardPropagationModel.hh
        G4INCLCascade.hh
        G4INCLCDPP.hh
        G4INCLDeltaProductionChannel.hh
        G4INCLCoulombNonRelativistic.hh
        G4INCLCrossSectionsMultiPions.hh
        G4INCLNuclearPotentialConstant.hh
        G4INCLNDFModifiedHarmonicOscillator.hh
        G4INCLPiNToDeltaChannel.hh
        G4INCLPauliStrict.hh
        G4INCLNuclearPotentialIsospin.hh
        G4INCLNNToMultiPionsChannel.hh
        G4INCLICoulomb.hh
        G4INCLPhaseSpaceGenerator.hh
        G4INCLCrossSections.hh
        G4INCLRecombinationChannel.hh
        G4INCLPauliStrictStandard.hh
        G4INCLNuclearDensity.hh
        G4INCLParticleSampler.hh
        G4INCLICrossSections.hh
        G4INCLConstantRandom.hh
        G4INCLNDFParis.hh
        G4INCLClusteringModelIntercomparison.hh

    SOURCES
        G4INCLPhaseSpaceKopylov.cc
        G4INCLPauliStandard.cc
        G4INCLParticleSampler.cc
        G4INCLNuclearDensityFactory.cc
        G4INCLPhaseSpaceRauboldLynch.cc
        G4INCLProjectileRemnant.cc
        G4INCLNuclearPotentialConstant.cc
        G4INCLCoulombDistortion.cc
        G4INCLStore.cc
        G4INCLReflectionChannel.cc
        G4INCLCoulombNone.cc
        G4INCLCDPP.cc
        G4INCLInteractionAvatar.cc
        G4INCLCascade.cc
        G4INCLINuclearPotential.cc
        G4INCLPiNElasticChannel.cc
        G4INCLParticleEntryChannel.cc
        G4INCLPauliBlocking.cc
        G4INCLDeuteronDensity.cc
        G4INCLCrossSectionsINCL46.cc
        G4INCLSurfaceAvatar.cc
        G4INCLDeltaDecayChannel.cc
        G4INCLNuclearPotentialIsospin.cc
        G4INCLNuclearPotentialEnergyIsospin.cc
        G4INCLRecombinationChannel.cc
        G4INCLParticleEntryAvatar.cc
        G4INCLDeltaProductionChannel.cc
        G4INCLNNToMultiPionsChannel.cc
        G4INCLDeJongSpin.cc
        G4INCLStandardPropagationModel.cc
        G4INCLElasticChannel.cc
        G4INCLPhaseSpaceGenerator.cc
        G4INCLDecayAvatar.cc
        G4INCLNuclearDensity.cc
        G4INCLClusterDecay.cc
        G4INCLAvatarDumpAction.cc
        G4INCLCrossSectionsMultiPions.cc
        G4INCLPiNToDeltaChannel.cc
        G4INCLCoulombNonRelativistic.cc
        G4INCLTransmissionChannel.cc
        G4INCLPauliStrict.cc
        G4INCLKinematicsUtils.cc
        G4INCLPauliGlobal.cc
        G4INCLPauliStrictStandard.cc
        G4INCLCrossSections.cc
        G4INCLPiNToMultiPionsChannel.cc
        G4INCLCascadeAction.cc
        G4INCLNucleus.cc
        G4INCLClusteringModelIntercomparison.cc
        G4INCLClustering.cc
        G4INCLCluster.cc
        G4INCLBinaryCollisionAvatar.cc
        G4INCLNuclearPotentialEnergyIsospinSmooth.cc

    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
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
        G4intercoms
        G4had_preequ_exciton
        G4hadronic_deex_evaporation
        G4hadronic_deex_fermi_breakup
        G4hadronic_deex_handler
        G4hadronic_deex_management
        G4hadronic_deex_multifragmentation
        G4hadronic_deex_photon_evaporation
        G4hadronic_deex_util
        G4hadronic_inclxx_utils

    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4materials
        G4particles
        G4track
        G4intercoms

LINK_LIBRARIES)# List any source specific properties here