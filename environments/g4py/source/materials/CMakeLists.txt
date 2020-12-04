# - build library

# library
geant4_add_pymodule(pyG4materials
  pyG4AtomicShells.cc
  pyG4Element.cc
  pyG4ElementTable.cc
  pyG4ElementVector.cc
  pyG4Isotope.cc
  pyG4Material.cc
  pyG4MaterialTable.cc
  pyG4NistManager.cc
  pymodG4materials.cc
)

target_link_libraries(pyG4materials PRIVATE G4materials)
install(TARGETS pyG4materials DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/Geant4")
