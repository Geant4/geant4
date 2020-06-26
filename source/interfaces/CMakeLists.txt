#------------------------------------------------------------------------------
# Module : G4interfaces
# Package: Geant4.src.G4interfaces
#------------------------------------------------------------------------------

geant4_global_library_target(NAME G4interfaces
  COMPONENTS
    GAG/sources.cmake
    basic/sources.cmake
    common/sources.cmake
)

if(GEANT4_USE_QT)
  foreach(__target G4interfaces G4interfaces-static)
    if(TARGET ${__target})
      set_target_properties(${__target} PROPERTIES AUTOMOC 1)
    endif()
  endforeach()
endif()