# library
geant4_add_pymodule(pyG4gdml
  pyG4GDMLParser.cc
  pymodG4g4gdml.cc
)

target_link_libraries(pyG4gdml PRIVATE G4persistency)
install(TARGETS pyG4gdml DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/Geant4")
