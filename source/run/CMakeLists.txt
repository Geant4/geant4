#------------------------------------------------------------------------------
# Module : G4run
# Package: Geant4.src.G4run
#------------------------------------------------------------------------------

# Add allocation export symbol for the run category
add_definitions(-DG4RUN_ALLOC_EXPORT)

geant4_global_library_target(COMPONENTS sources.cmake)
