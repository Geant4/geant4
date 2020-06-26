#------------------------------------------------------------------------------
# Module : G4particles
# Package: Geant4.src.G4particles
#------------------------------------------------------------------------------

# Add allocation export symbol for the particles category
add_definitions(-DG4PARTICLES_ALLOC_EXPORT)

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
    utils/sources.cmake
)

