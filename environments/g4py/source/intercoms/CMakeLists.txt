# - build library

# library
geant4_add_pymodule(pyG4intercoms
  pyG4UIcommand.cc
  pyG4UIcommandTree.cc
  pyG4UImanager.cc
  pyG4UIparameter.cc
  pymodG4intercoms.cc
)

target_link_libraries(pyG4intercoms PRIVATE G4intercoms)
install(TARGETS pyG4intercoms DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/Geant4")
