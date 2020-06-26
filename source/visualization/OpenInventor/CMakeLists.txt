#------------------------------------------------------------------------------
# Module : G4OpenInventor
# Package: Geant4.src.G4visualization.G4OpenInventor
#------------------------------------------------------------------------------

geant4_global_library_target(COMPONENTS sources.cmake)

# Fix for Inventor not yet supporting usage requirements
if(TARGET G4OpenInventor)
  target_include_directories(G4OpenInventor PUBLIC ${G4VIS_MODULE_OPENINVENTOR_INCLUDE_DIRS})
endif()

if(TARGET G4OpenInventor-static)
  target_include_directories(G4OpenInventor-static PUBLIC ${G4VIS_MODULE_OPENINVENTOR_INCLUDE_DIRS})
endif()

if(GEANT4_USE_QT)
  foreach(__target G4OpenInventor G4OpenInventor-static)
    if(TARGET ${__target})
      set_target_properties(${__target} PROPERTIES AUTOMOC 1)
    endif()
  endforeach()
endif()