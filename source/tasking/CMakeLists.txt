#-----------------------------------------------------------------------
# Module : G4tasking
# Package: Geant4.src.G4tasking
#-----------------------------------------------------------------------

# Add allocation export symbol for the tasking category
add_definitions(-DG4TASKING_ALLOC_EXPORT)

geant4_global_library_target(NAME G4tasking
  COMPONENTS
    sources.cmake
  )

if(TARGET G4tasking)
  target_include_directories(G4tasking PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source/externals/ptl>
  )
endif()

if(TARGET G4tasking-static)
  target_include_directories(G4tasking-static PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source/externals/ptl>
  )
endif()

if(GEANT4_USE_TBB)
    foreach(_LIB G4tasking G4tasking-static)
        if(TARGET ${_LIB})
            target_include_directories(${_LIB} PUBLIC ${TBB_INCLUDE_DIRS})
            target_link_libraries(${_LIB} ${TBB_LIBRARIES})
        endif()
    endforeach()
endif()
