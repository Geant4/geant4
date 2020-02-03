#-----------------------------------------------------------------------
# CMakeLists.txt
# Module : G4global
# Package: Geant4.src.G4global
#
# Top level CMakeLists.txt for a Geant4 Global Library
# or collection of sublibraries
#
# Generated on : 24/9/2010
#
#
#-----------------------------------------------------------------------

# Add allocation export symbol for the global category
add_definitions(-DG4GLOB_ALLOC_EXPORT)

add_subdirectory(HEPGeometry)
add_subdirectory(HEPNumerics)
add_subdirectory(HEPRandom)
add_subdirectory(management)

# Geant4 global library
if(NOT GEANT4_BUILD_GRANULAR_LIBS)
  include(Geant4MacroLibraryTargets)
  GEANT4_GLOBAL_LIBRARY_TARGET(NAME G4global
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
endif()
