# - G4phys_ctor_stopping module build definition

# Define the Geant4 Module.
geant4_add_module(G4phys_ctor_stopping
  PUBLIC_HEADERS
    G4StoppingPhysics.hh
    G4StoppingPhysicsFritiofWithBinaryCascade.hh
  SOURCES
    G4StoppingPhysics.cc
    G4StoppingPhysicsFritiofWithBinaryCascade.cc)

geant4_module_link_libraries(G4phys_ctor_stopping
  PUBLIC
    G4globman
    G4run
  PRIVATE
    G4baryons
    G4hadronic_stop
    G4leptons
    G4mesons
    G4partman
    G4phys_builders
    G4phys_ctor_factory
    G4procman)
