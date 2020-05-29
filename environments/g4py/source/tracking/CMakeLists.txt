# - build library

# library
geant4_add_pymodule(pyG4tracking
  pyG4TrackingManager.cc
  pyG4UserSteppingAction.cc
  pyG4UserTrackingAction.cc
  pymodG4tracking.cc
)

target_link_libraries(pyG4tracking PRIVATE G4tracking)
install(TARGETS pyG4tracking DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/Geant4")
