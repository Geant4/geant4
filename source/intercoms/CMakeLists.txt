#------------------------------------------------------------------------------
# Module : G4intercoms
# Package: Geant4.src.G4intercoms
#------------------------------------------------------------------------------

# Add allocation export symbol for the intercoms category
add_definitions(-DG4ICOMS_ALLOC_EXPORT)

geant4_global_library_target(COMPONENTS sources.cmake)
