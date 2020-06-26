#------------------------------------------------------------------------------
# Module : G4OpenGL
# Package: Geant4.src.G4visualization.G4OpenGL
#------------------------------------------------------------------------------

geant4_global_library_target(COMPONENTS sources.cmake)

if(GEANT4_USE_QT)
  foreach(__target G4OpenGL G4OpenGL-static)
    if(TARGET ${__target})
      set_target_properties(${__target} PROPERTIES AUTOMOC 1)
    endif()
  endforeach()
endif()