# - G4phys_ctor_ions module build definition

# Define the Geant4 Module.
geant4_add_module(G4phys_ctor_ions
  PUBLIC_HEADERS
    G4IonBinaryCascadePhysics.hh
    G4IonINCLXXPhysics.hh
    G4IonPhysics.hh
    G4IonPhysicsPHP.hh
    G4IonPhysicsXS.hh
    G4IonQMDPhysics.hh
    G4LightIonQMDPhysics.hh
  SOURCES
    G4IonBinaryCascadePhysics.cc
    G4IonINCLXXPhysics.cc
    G4IonPhysics.cc
    G4IonPhysicsPHP.cc
    G4IonPhysicsXS.cc
    G4IonQMDPhysics.cc
    G4LightIonQMDPhysics.cc)

geant4_module_link_libraries(G4phys_ctor_ions
  PUBLIC
    G4globman
    G4run
  PRIVATE
    G4had_par_hp
    G4had_preequ_exciton
    G4hadronic_binary
    G4hadronic_deex_handler
    G4hadronic_deex_management
    G4hadronic_inclxx_interface
    G4hadronic_mgt
    G4hadronic_proc
    G4hadronic_qmd
    G4hadronic_util
    G4hadronic_xsect
    G4ions
    G4partman
    G4phys_builders
    G4phys_ctor_factory
    G4procman)
