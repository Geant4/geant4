# - G4track category build

# Add (private) allocation export symbol
add_definitions(-DG4TRACK_ALLOC_EXPORT)

geant4_global_library_target(COMPONENTS sources.cmake)
