# - build library

# library
geant4_add_pymodule(pyG4physicslists
  pyPhysicsLists.cc
  pyPhysListFactory.cc
  pymodG4physicslists.cc
)
target_link_libraries(pyG4physicslists PRIVATE G4physicslists)
install(TARGETS pyG4physicslists DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/Geant4")
