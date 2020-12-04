# - build library

# library
geant4_add_pymodule(pyG4processes
  pyG4CrossSectionHandler.cc
  pyG4EmCalculator.cc
  pyG4ProcVector.cc
  pyG4ProcessManager.cc
  pyG4ProcessTable.cc
  pyG4ProcessType.cc
  pyG4ProductionCutsTable.cc
  pyG4VCrossSectionHandler.cc
  pyG4VProcess.cc
  pymodG4processes.cc
)

target_link_libraries(pyG4processes PRIVATE G4processes)
install(TARGETS pyG4processes DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/Geant4")
