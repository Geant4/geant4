# - G4processes category build

# Add (private) allocation export symbol for the EM category
add_definitions(-DG4EM_ALLOC_EXPORT)

# Add (private) allocation export symbol for the hadronic category
add_definitions(-DG4HADRONIC_ALLOC_EXPORT)

geant4_global_library_target(NAME G4processes
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
    hadronic/models/gamma_nuclear/sources.cmake
    hadronic/models/lend/sources.cmake
    hadronic/models/lepto_nuclear/sources.cmake
    hadronic/models/particle_hp/sources.cmake
    hadronic/models/parton_string/diffraction/sources.cmake
    hadronic/models/parton_string/hadronization/sources.cmake
    hadronic/models/parton_string/management/sources.cmake
    hadronic/models/parton_string/qgsm/sources.cmake
    hadronic/models/pre_equilibrium/exciton_model/sources.cmake
    hadronic/models/qmd/sources.cmake
    hadronic/models/quasi_elastic/sources.cmake
    hadronic/models/radioactive_decay/sources.cmake
    hadronic/models/theo_high_energy/sources.cmake
    hadronic/processes/sources.cmake
    hadronic/stopping/sources.cmake
    hadronic/util/sources.cmake
    management/sources.cmake
    optical/sources.cmake
    solidstate/phonon/sources.cmake
    solidstate/channeling/sources.cmake
    parameterisation/sources.cmake
    scoring/sources.cmake
    transportation/sources.cmake)

