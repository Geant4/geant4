# - Build for G4tools as an "external"
# It can be built as a standard module, but we don't declare the headers
# so that their nested structure can be retained on install
geant4_add_module(G4tools SOURCES dummy.cc)

geant4_add_category(G4tools MODULES G4tools)

#-----------------------------------------------------------------------
# Custom install of (nested) headers and license
install(DIRECTORY include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
  COMPONENT Development
  )

install(FILES tools.license DESTINATION ${CMAKE_INSTALL_DATADIR})

#-----------------------------------------------------------------------
# Install fonts
#
if(GEANT4_USE_FREETYPE)
  install(DIRECTORY fonts/ DESTINATION ${CMAKE_INSTALL_DATADIR}/fonts)
endif()
