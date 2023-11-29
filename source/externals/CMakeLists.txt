# - Build of external packages, if required

if(NOT GEANT4_USE_SYSTEM_CLHEP)
  add_subdirectory(clhep)
endif()

if(NOT GEANT4_USE_SYSTEM_EXPAT)
  add_subdirectory(expat)
endif()

if(NOT GEANT4_USE_SYSTEM_ZLIB)
  add_subdirectory(zlib)
endif()

add_subdirectory(g4tools)

if(NOT GEANT4_USE_SYSTEM_PTL)
  # We reuse PTL as a true subproject, configure settings before recursion
  set(PTL_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
  set(PTL_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}/PTL)
  set(PTL_INSTALL_HEADERS ON)
  set(PTL_INSTALL_CONFIG ON)
  set(PTL_USE_TBB ${GEANT4_USE_TBB})
  set(PTL_USE_LOCKS ${GEANT4_USE_PTL_LOCKS})

  if(GEANT4_BUILD_SANITIZER AND NOT (GEANT4_BUILD_SANITIZER STREQUAL "none"))
    set(PTL_USE_SANITIZER ON)
    set(PTL_SANITIZER_TYPE ${GEANT4_BUILD_SANITIZER})
  endif()

  add_subdirectory(ptl)

  # Ensure include path is available from build dir in scripts
  set_property(GLOBAL APPEND PROPERTY GEANT4_BUILDTREE_INCLUDE_DIRS 
    ${CMAKE_CURRENT_SOURCE_DIR}/ptl/include)

  # Adapt targets to G4 naming scheme, whether
  if(TARGET ptl-shared)
    add_library(G4ptl INTERFACE)
    target_link_libraries(G4ptl INTERFACE PTL::ptl-shared)
    install(TARGETS G4ptl
      EXPORT Geant4LibraryDepends
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
      INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
    add_library(Geant4::G4ptl ALIAS G4ptl)
  endif()

  if(TARGET ptl-static)
    add_library(G4ptl-static INTERFACE)
    target_link_libraries(G4ptl-static INTERFACE PTL::ptl-static)
    install(TARGETS G4ptl-static
      EXPORT Geant4LibraryDepends
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
      INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
    add_library(Geant4::G4ptl-static ALIAS G4ptl-static)
  endif()

  geant4_add_external_category(G4ptl)
endif()
