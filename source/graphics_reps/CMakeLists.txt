#------------------------------------------------------------------------------
# Module : G4graphics_reps
# Package: Geant4.src.G4graphics_reps
#------------------------------------------------------------------------------

# Add allocation export symbol for the graphics_reps category
add_definitions(-DG4GREPS_ALLOC_EXPORT)

geant4_global_library_target(COMPONENTS sources.cmake)
