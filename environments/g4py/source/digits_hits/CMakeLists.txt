# - build library

# library
geant4_add_pymodule(pyG4digits_hits
  pyG4VSensitiveDetector.cc
  pymodG4digits_hits.cc
)

target_link_libraries(pyG4digits_hits PRIVATE G4digits_hits)
install(TARGETS pyG4digits_hits DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/Geant4")
