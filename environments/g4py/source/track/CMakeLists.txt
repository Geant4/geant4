# - build library

# library
geant4_add_pymodule(pyG4track
  pyG4Step.cc
  pyG4StepPoint.cc
  pyG4StepStatus.cc
  pyG4Track.cc
  pyG4TrackStatus.cc
  pymodG4track.cc
)

target_link_libraries(pyG4track PRIVATE G4track G4processes)
install(TARGETS pyG4track DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/Geant4")
