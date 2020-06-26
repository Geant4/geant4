#------------------------------------------------------------------------------
# Module : G4digits_hits
# Package: Geant4.src.G4digits_hits
#------------------------------------------------------------------------------

# Add allocation export symbol for the digits_hits category
add_definitions(-DG4DIGI_ALLOC_EXPORT)

geant4_global_library_target(NAME G4digits_hits
  COMPONENTS
    detector/sources.cmake
    digits/sources.cmake
    hits/sources.cmake
    scorer/sources.cmake
    utils/sources.cmake
)
