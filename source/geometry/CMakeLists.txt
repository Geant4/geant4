# - G4geometry category builds

# Add (private) allocation export symbol
add_definitions(-DG4GEOM_ALLOC_EXPORT)

# Compose library from subcategories
geant4_global_library_target(NAME G4geometry
  COMPONENTS
    biasing/sources.cmake
    divisions/sources.cmake
    magneticfield/sources.cmake
    management/sources.cmake
    navigation/sources.cmake
    solids/Boolean/sources.cmake
    solids/CSG/sources.cmake
    solids/specific/sources.cmake
    volumes/sources.cmake)
