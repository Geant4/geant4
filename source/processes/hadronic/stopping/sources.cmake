#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_stop
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_stop
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake,v 1.1 2010-09-29 19:08:23 bmorgan Exp $
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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/evaporation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/fermi_breakup/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/fission/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/gem_evaporation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/handler/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/multifragmentation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/photon_evaporation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/de_excitation/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/high_energy/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/pre_equilibrium/exciton_model/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/models/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/hadronic/util/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hadronic_stop
    HEADERS
        G4AntiNeutronAnnihilationAtRest.hh
        G4AntiProtonAnnihilationAtRest.hh
        G4DistributionGenerator.hh
        G4KaonMinusAbsorption.hh
        G4KaonMinusAbsorptionAtRest.hh
        G4MuMinusCaptureCascade.hh
        G4MuonMinusCaptureAtRest.hh
        G4NeutronCaptureAtRest.hh
        G4PiMinusAbsorptionAtRest.hh
        G4PiMinusStopAbsorption.hh
        G4PiMinusStopAl.hh
        G4PiMinusStopC.hh
        G4PiMinusStopCo.hh
        G4PiMinusStopCu.hh
        G4PiMinusStopLi.hh
        G4PiMinusStopMaterial.hh
        G4PiMinusStopN.hh
        G4PiMinusStopO.hh
        G4PiMinusStopPb.hh
        G4PiMinusStopTa.hh
        G4PionMinusAbsorptionAtRest.hh
        G4StopDeexcitation.hh
        G4StopDeexcitationAlgorithm.hh
        G4StopDummyDeexcitation.hh
        G4StopElementSelector.hh
        G4StopTheoDeexcitation.hh
    SOURCES
        G4AntiNeutronAnnihilationAtRest.cc
        G4AntiProtonAnnihilationAtRest.cc
        G4DistributionGenerator.cc
        G4KaonMinusAbsorption.cc
        G4KaonMinusAbsorptionAtRest.cc
        G4MuMinusCaptureCascade.cc
        G4MuonMinusCaptureAtRest.cc
        G4NeutronCaptureAtRest.cc
        G4PiMinusAbsorptionAtRest.cc
        G4PiMinusStopAbsorption.cc
        G4PiMinusStopAl.cc
        G4PiMinusStopC.cc
        G4PiMinusStopCo.cc
        G4PiMinusStopCu.cc
        G4PiMinusStopLi.cc
        G4PiMinusStopMaterial.cc
        G4PiMinusStopN.cc
        G4PiMinusStopO.cc
        G4PiMinusStopPb.cc
        G4PiMinusStopTa.cc
        G4PionMinusAbsorptionAtRest.cc
        G4StopDeexcitation.cc
        G4StopDummyDeexcitation.cc
        G4StopElementSelector.cc
        G4StopTheoDeexcitation.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
        G4had_mod_man
        G4had_mod_util
        G4had_preequ_exciton
        G4hadronic_HE
        G4hadronic_deex_evaporation
        G4hadronic_deex_fermi_breakup
        G4hadronic_deex_fission
        G4hadronic_deex_gem_evaporation
        G4hadronic_deex_handler
        G4hadronic_deex_management
        G4hadronic_deex_multifragmentation
        G4hadronic_deex_photon_evaporation
        G4hadronic_deex_util
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

