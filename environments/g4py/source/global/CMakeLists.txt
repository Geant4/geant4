# - Module
geant4_add_pymodule(pyG4global
  G4PyCoutDestination.cc
  pyG4ApplicationState.cc
  pyG4Exception.cc
  pyG4ExceptionHandler.cc
  pyG4ExceptionSeverity.cc
  pyG4RandomDirection.cc
  pyG4RotationMatrix.cc
  pyG4StateManager.cc
  pyG4String.cc
  pyG4ThreeVector.cc
  pyG4Timer.cc
  pyG4Transform3D.cc
  pyG4TwoVector.cc
  pyG4UnitsTable.cc
  pyG4UserLimits.cc
  pyG4Version.cc
  pyRandomEngines.cc
  pyRandomize.cc
  pygeomdefs.cc
  pyglobals.cc
  pymodG4global.cc
)

target_link_libraries(pyG4global PRIVATE G4global G4track)
install(TARGETS pyG4global DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/Geant4")
