# - Build for G4tools as an external
# As a pure interface, we simply declare and install it as such.
# It needs no coupling to the internal/external category system because it is header-only.
add_library(G4tools INTERFACE)
add_library(Geant4::G4tools ALIAS G4tools)
target_compile_features(G4tools INTERFACE ${GEANT4_TARGET_COMPILE_FEATURES})
target_include_directories(G4tools 
  INTERFACE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}> )

# Custom install of target, (nested) headers, license and fonts
install(TARGETS G4tools EXPORT Geant4LibraryDepends)  
install(DIRECTORY include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
  COMPONENT Development )
install(FILES tools.license DESTINATION ${CMAKE_INSTALL_DATADIR})
install(DIRECTORY fonts/ DESTINATION ${CMAKE_INSTALL_DATADIR}/fonts)
