# - G4intercoms category build(s)

# Add allocation export symbol for the intercoms category
add_definitions(-DG4ICOMS_ALLOC_EXPORT)

geant4_global_library_target(NAME G4intercoms COMPONENTS sources.cmake)
