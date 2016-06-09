#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_bert_cascade
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4had_hetcpp.G4hadronic_bert_cascade
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 30/9/2010
#
# $Id: sources.cmake,v 1.4 2010-09-30 12:02:28 bmorgan Exp $
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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/cascade/cascade/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/cascade/utils/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/evaporation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/fermi_breakup/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/handler/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/multifragmentation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/photon_evaporation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/pre_equilibrium/exciton_model/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/processes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_inclxx_physics
    HEADERS
		   G4INCLAvatarAction.hh
			 G4INCLBinaryCollisionAvatar.hh
			 G4INCLCDPP.hh
			 G4INCLClusterDecay.hh
			 G4INCLClustering.hh
			 G4INCLClusteringModelIntercomparison.hh
			 G4INCLClusteringModelNone.hh
			 G4INCLClusterUtils.hh
			 G4INCLConstantRandom.hh
			 G4INCLCoulombDistortion.hh
			 G4INCLCoulombNone.hh
			 G4INCLCoulombNonRelativistic.hh
			 G4INCLCrossSections.hh
			 G4INCLDecayAvatar.hh
			 G4INCLDeltaDecayChannel.hh
			 G4INCLDeltaProductionChannel.hh
			 G4INCLElasticChannel.hh
			 G4INCLEventAction.hh
			 G4INCLFPEDebug.hh
			 G4INCLGlobalInfo.hh
			 G4INCLCascade.hh
			 G4INCLIClusteringModel.hh
			 G4INCLICoulomb.hh
			 G4INCLIFunction.hh
			 G4INCLInteractionAvatar.hh
			 G4INCLINuclearPotential.hh
			 G4INCLIPauli.hh
			 G4INCLIPropagationModel.hh
			 G4INCLKinematicsUtils.hh
			 G4INCLNuclearDensityFactory.hh
			 G4INCLNuclearDensity.hh
			 G4INCLNuclearPotentialConstant.hh
			 G4INCLNuclearPotentialEnergyIsospin.hh
			 G4INCLNuclearPotentialIsospin.hh
			 G4INCLNucleus.hh
			 G4INCLPauliBlocking.hh
			 G4INCLPauliGlobal.hh
			 G4INCLPauliStandard.hh
			 G4INCLPauliStrict.hh
			 G4INCLPauliStrictStandard.hh
			 G4INCLPionNucleonChannel.hh
			 G4INCLPropagationAction.hh
			 G4INCLRecombinationChannel.hh
			 G4INCLReflectionChannel.hh
			 G4INCLStandardPropagationModel.hh
			 G4INCLSurfaceAvatar.hh
			 G4INCLTransmissionChannel.hh
    SOURCES
		   G4INCLAvatarAction.cc
			 G4INCLBinaryCollisionAvatar.cc
			 G4INCLCascade.cc
			 G4INCLCDPP.cc
			 G4INCLClusterDecay.cc
			 G4INCLClustering.cc
			 G4INCLClusteringModelIntercomparison.cc
			 G4INCLClusterUtils.cc
			 G4INCLCoulombDistortion.cc
			 G4INCLCoulombNone.cc
			 G4INCLCoulombNonRelativistic.cc
			 G4INCLCrossSections.cc
			 G4INCLDecayAvatar.cc
			 G4INCLDeltaDecayChannel.cc
			 G4INCLDeltaProductionChannel.cc
			 G4INCLElasticChannel.cc
			 G4INCLEventAction.cc
			 G4INCLIClusteringModel.cc
			 G4INCLICoulomb.cc
			 G4INCLInteractionAvatar.cc
			 G4INCLINuclearPotential.cc
			 G4INCLIPropagationModel.cc
			 G4INCLKinematicsUtils.cc
			 G4INCLNuclearDensity.cc
			 G4INCLNuclearDensityFactory.cc
			 G4INCLNuclearPotentialConstant.cc
			 G4INCLNuclearPotentialEnergyIsospin.cc
			 G4INCLNuclearPotentialIsospin.cc
			 G4INCLNucleus.cc
			 G4INCLPauliBlocking.cc
			 G4INCLPauliGlobal.cc
			 G4INCLPauliStandard.cc
			 G4INCLPauliStrict.cc
			 G4INCLPauliStrictStandard.cc
			 G4INCLPionNucleonChannel.cc
			 G4INCLPropagationAction.cc
			 G4INCLRecombinationChannel.cc
			 G4INCLReflectionChannel.cc
			 G4INCLStandardPropagationModel.cc
			 G4INCLSurfaceAvatar.cc
			 G4INCLTransmissionChannel.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4had_mod_man
        G4had_mod_util
        G4had_preequ_exciton
        G4hadronic_inclxx_utils
        G4hadronic_deex_evaporation
        G4hadronic_deex_fermi_breakup
        G4hadronic_deex_handler
        G4hadronic_deex_management
        G4hadronic_deex_multifragmentation
        G4hadronic_deex_photon_evaporation
        G4hadronic_deex_util
        G4hadronic_hetcpp_utils
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

