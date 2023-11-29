# - G4particles category build
add_definitions(-DG4PARTICLES_ALLOC_EXPORT)

# Compose library from subcategories
foreach(modfile
    adjoint/sources.cmake
    bosons/sources.cmake
    hadrons/barions/sources.cmake
    hadrons/ions/sources.cmake
    hadrons/mesons/sources.cmake
    leptons/sources.cmake
    management/sources.cmake
    shortlived/sources.cmake
    utils/sources.cmake)
  include(${modfile})
endforeach()

geant4_add_category(G4particles
  MODULES
    G4partadj
    G4bosons
    G4baryons
    G4ions
    G4mesons
    G4leptons
    G4partman
    G4shortlived
    G4partutils )
