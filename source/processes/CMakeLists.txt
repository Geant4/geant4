#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4processes
# Package: Geant4.src.G4processes
#
# Top level CMakeLists.txt for a Geant4 Global Library
# or collection of sublibraries
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 103639 2017-04-20 09:03:42Z gcosmo $
#
# 20130225  M. Kelsey -- Remove cascade/evaporation, cascade/utils
#------------------------------------------------------------------------------

# Add allocation export symbol for the EM category
add_definitions(-DG4EM_ALLOC_EXPORT)

# Add allocation export symbol for the hadronic category
add_definitions(-DG4HADRONIC_ALLOC_EXPORT)

add_subdirectory(biasing)
add_subdirectory(cuts)
add_subdirectory(decay)
add_subdirectory(electromagnetic)
add_subdirectory(hadronic)
add_subdirectory(management)
add_subdirectory(optical)
add_subdirectory(solidstate)
add_subdirectory(parameterisation)
add_subdirectory(scoring)
add_subdirectory(transportation)

if(NOT GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GLOBAL_LIBRARY_TARGET(NAME G4processes
        COMPONENTS
            biasing/management/sources.cmake
            biasing/generic/sources.cmake
            biasing/importance/sources.cmake
            cuts/sources.cmake
            decay/sources.cmake
            electromagnetic/adjoint/sources.cmake
            electromagnetic/dna/processes/sources.cmake
            electromagnetic/dna/models/sources.cmake
            electromagnetic/dna/utils/sources.cmake
            electromagnetic/dna/management/sources.cmake
            electromagnetic/dna/molecules/management/sources.cmake
            electromagnetic/dna/molecules/types/sources.cmake
            electromagnetic/highenergy/sources.cmake
            electromagnetic/lowenergy/sources.cmake
            electromagnetic/muons/sources.cmake
            electromagnetic/pii/sources.cmake
            electromagnetic/polarisation/sources.cmake
            electromagnetic/standard/sources.cmake
            electromagnetic/utils/sources.cmake
            electromagnetic/xrays/sources.cmake
            hadronic/cross_sections/sources.cmake
            hadronic/management/sources.cmake
            hadronic/models/abla/sources.cmake
            hadronic/models/abrasion/sources.cmake
            hadronic/models/binary_cascade/sources.cmake
            hadronic/models/cascade/cascade/sources.cmake
            hadronic/models/coherent_elastic/sources.cmake
            hadronic/models/de_excitation/ablation/sources.cmake
            hadronic/models/de_excitation/evaporation/sources.cmake
            hadronic/models/de_excitation/fermi_breakup/sources.cmake
            hadronic/models/de_excitation/fission/sources.cmake
            hadronic/models/de_excitation/gem_evaporation/sources.cmake
            hadronic/models/de_excitation/handler/sources.cmake
            hadronic/models/de_excitation/management/sources.cmake
            hadronic/models/de_excitation/multifragmentation/sources.cmake
            hadronic/models/de_excitation/photon_evaporation/sources.cmake
            hadronic/models/de_excitation/util/sources.cmake
            hadronic/models/em_dissociation/sources.cmake
            hadronic/models/fission/sources.cmake
            hadronic/models/im_r_matrix/sources.cmake
            hadronic/models/inclxx/utils/sources.cmake
            hadronic/models/inclxx/incl_physics/sources.cmake
            hadronic/models/inclxx/interface/sources.cmake
            hadronic/models/lend/sources.cmake
            hadronic/models/lepto_nuclear/sources.cmake
            hadronic/models/management/sources.cmake
            hadronic/models/particle_hp/sources.cmake
            hadronic/models/parton_string/diffraction/sources.cmake
            hadronic/models/parton_string/hadronization/sources.cmake
            hadronic/models/parton_string/management/sources.cmake
            hadronic/models/parton_string/qgsm/sources.cmake
            hadronic/models/pre_equilibrium/exciton_model/sources.cmake
            hadronic/models/qmd/sources.cmake
            hadronic/models/radioactive_decay/sources.cmake
            hadronic/models/quasi_elastic/sources.cmake
            hadronic/models/rpg/sources.cmake
            hadronic/models/theo_high_energy/sources.cmake
            hadronic/models/util/sources.cmake
            hadronic/processes/sources.cmake
            hadronic/stopping/sources.cmake
            hadronic/util/sources.cmake
            management/sources.cmake
            optical/sources.cmake
            solidstate/phonon/sources.cmake
            solidstate/channeling/sources.cmake
            parameterisation/sources.cmake
            scoring/sources.cmake
            transportation/sources.cmake
)
endif()

