# - build library

# library
geant4_add_pymodule(pyG4graphics_reps
  pyG4Colour.cc
  pyG4VisAttributes.cc
  pymodG4graphics_reps.cc
)

target_link_libraries(pyG4graphics_reps PRIVATE G4graphics_reps)
install(TARGETS pyG4graphics_reps DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/Geant4")
