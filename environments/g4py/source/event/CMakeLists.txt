# - build library

# library
geant4_add_pymodule(pyG4event
  pyG4ClassificationOfNewTrack.cc
  pyG4Event.cc
  pyG4EventManager.cc
  pyG4ParticleGun.cc
  pyG4StackManager.cc
  pyG4UserEventAction.cc
  pyG4UserStackingAction.cc
  pymodG4event.cc
)

target_link_libraries(pyG4event PRIVATE G4event)
install(TARGETS pyG4event DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/Geant4")
