#------------------------------------------------------------------------------
# Module : G4tracking
# Package: Geant4.src.G4tracking
#------------------------------------------------------------------------------

# Add allocation export symbol for the tracking category
add_definitions(-DG4TRACKING_ALLOC_EXPORT)

geant4_global_library_target(COMPONENTS sources.cmake)
