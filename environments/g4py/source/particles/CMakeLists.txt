# - build library

# library
geant4_add_pymodule(pyG4particles
  pyG4DecayTable.cc
  pyG4DynamicParticle.cc
  pyG4ParticleDefinition.cc
  pyG4ParticleList.cc
  pyG4ParticleTable.cc
  pyG4PrimaryParticle.cc
  pyG4PrimaryVertex.cc
  pymodG4particles.cc
)

target_link_libraries(pyG4particles PRIVATE G4particles G4processes)
install(TARGETS pyG4particles DESTINATION "${CMAKE_INSTALL_PYTHONDIR}/Geant4")
