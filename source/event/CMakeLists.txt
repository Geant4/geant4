#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4event
# Package: Geant4.src.G4event
#------------------------------------------------------------------------------

# Add allocation export symbol for the event category
add_definitions(-DG4EVENT_ALLOC_EXPORT)

geant4_global_library_target(COMPONENTS sources.cmake)

# Enable smart stack
if(GEANT4_USE_SMARTSTACK)
  if(TARGET G4event)
    target_compile_definitions(G4event PUBLIC G4_USESMARTSTACK)
  endif()

  if(TARGET G4event-static)
    target_compile_definitions(G4event-static PUBLIC G4_USESMARTSTACK)
  endif()
endif()