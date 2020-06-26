#------------------------------------------------------------------------------
# Module : G4geometry
# Package: Geant4.src.G4geometry
#------------------------------------------------------------------------------

# Add allocation export symbol for the geometry category
add_definitions(-DG4GEOM_ALLOC_EXPORT)

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
    volumes/sources.cmake
  )

# Temp hack to get "G4GeomConfig.hh" generated header in include
# path for dependent targets
# TODO: Remove once migration to new module/library system in place.
#       This will automatically set include paths
if(TARGET G4geometry)
  target_include_directories(G4geometry PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
    )
endif()

if(TARGET G4geometry-static)
  target_include_directories(G4geometry-static PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  )
endif()

