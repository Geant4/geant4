# - G4event category build

# Add (private) allocation export symbol
add_definitions(-DG4EVENT_ALLOC_EXPORT)

geant4_global_library_target(COMPONENTS sources.cmake)
