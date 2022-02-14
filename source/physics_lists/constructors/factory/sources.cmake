# - G4phys_ctor_factory module build definition

# Define the Geant4 Module.
geant4_add_module(G4phys_ctor_factory
  PUBLIC_HEADERS
    G4RegisterPhysicsConstructors.icc
    G4PhysicsConstructorFactory.hh
    G4PhysicsConstructorRegistry.hh
  SOURCES
    G4PhysicsConstructorRegistry.cc)

geant4_module_link_libraries(G4phys_ctor_factory PUBLIC G4globman G4run)
