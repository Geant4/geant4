# - G4digits_hits category build

# Add (private) allocation export symbol
add_definitions(-DG4DIGI_ALLOC_EXPORT)

geant4_global_library_target(NAME G4digits_hits
  COMPONENTS
    detector/sources.cmake
    digits/sources.cmake
    hits/sources.cmake
    scorer/sources.cmake
    utils/sources.cmake)
