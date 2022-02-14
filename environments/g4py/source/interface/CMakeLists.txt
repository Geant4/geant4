# - build library

# library
geant4_add_pymodule(pyG4interfaces
  pyG4UIterminal.cc
  pymodG4interface.cc
)

target_link_libraries(pyG4interfaces PRIVATE G4interfaces)
install(TARGETS pyG4interfaces DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/Geant4")
