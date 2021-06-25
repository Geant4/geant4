# - G4run category build

# Add (private) allocation export symbol
add_definitions(-DG4RUN_ALLOC_EXPORT)

geant4_global_library_target(COMPONENTS sources.cmake)
