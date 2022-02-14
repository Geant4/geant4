# - G4global category build(s)

# Add (private) allocation export symbol
add_definitions(-DG4GLOB_ALLOC_EXPORT)

# Compose library from subcategories
geant4_global_library_target(NAME G4global
  COMPONENTS
    HEPGeometry/sources.cmake
    HEPNumerics/sources.cmake
    HEPRandom/sources.cmake
    management/sources.cmake
  )
