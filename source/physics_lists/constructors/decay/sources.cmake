# - G4phys_ctor_decay module build definition

# Define the Geant4 Module.
geant4_add_module(G4phys_ctor_decay
  PUBLIC_HEADERS
    G4DecayPhysics.hh
    G4MuonicAtomDecayPhysics.hh
    G4SpinDecayPhysics.hh
    G4RadioactiveDecayPhysics.hh
    G4UnknownDecayPhysics.hh
  SOURCES
    G4DecayPhysics.cc
    G4MuonicAtomDecayPhysics.cc
    G4SpinDecayPhysics.cc
    G4RadioactiveDecayPhysics.cc
    G4UnknownDecayPhysics.cc)

geant4_module_link_libraries(G4phys_ctor_decay
  PUBLIC
    G4decay
    G4globman
    G4run
  PRIVATE
    G4baryons
    G4bosons
    G4emlowenergy
    G4emutils
    G4hadronic_deex_management
    G4hadronic_radioactivedecay
    G4hadronic_stop
    G4ions
    G4leptons
    G4mesons
    G4partman
    G4phys_builders
    G4phys_ctor_em
    G4phys_ctor_factory
    G4physlist_util
    G4procman
    G4shortlived)
