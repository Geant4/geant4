# - G4graphics_reps category build(s)

# Add allocation export symbol for the graphics_reps category
add_definitions(-DG4GREPS_ALLOC_EXPORT)

geant4_global_library_target(NAME G4graphics_reps COMPONENTS sources.cmake)
