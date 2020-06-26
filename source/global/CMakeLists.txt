#-----------------------------------------------------------------------
# Module : G4global
# Package: Geant4.src.G4global
#-----------------------------------------------------------------------

# Add allocation export symbol for the global category
add_definitions(-DG4GLOB_ALLOC_EXPORT)

geant4_global_library_target(NAME G4global
  COMPONENTS
    HEPGeometry/sources.cmake
    HEPNumerics/sources.cmake
    HEPRandom/sources.cmake
    management/sources.cmake
  )

# Temp hack to get "G4GlobalConfig.hh" generated header in include
# path for dependent targets
# TODO: Remove once migration to new module/library system in place.
#       This will automatically set include paths
if(TARGET G4global)
  target_include_directories(G4global PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
    )
endif()

if(TARGET G4global-static)
  target_include_directories(G4global-static PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  )
endif()
