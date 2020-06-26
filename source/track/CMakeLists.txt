#------------------------------------------------------------------------------
# Module : G4track
# Package: Geant4.src.G4track
#------------------------------------------------------------------------------

# Add allocation export symbol for the track category
add_definitions(-DG4TRACK_ALLOC_EXPORT)

geant4_global_library_target(COMPONENTS sources.cmake)
