# - G4particles category build

# Add (private) allocation export symbol
add_definitions(-DG4PARTICLES_ALLOC_EXPORT)

# Compose library from subcategories
geant4_global_library_target(NAME G4particles
  COMPONENTS
    adjoint/sources.cmake
    bosons/sources.cmake
    hadrons/barions/sources.cmake
    hadrons/ions/sources.cmake
    hadrons/mesons/sources.cmake
    leptons/sources.cmake
    management/sources.cmake
    shortlived/sources.cmake
    utils/sources.cmake)
