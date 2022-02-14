# - G4hadronic_inclxx_interface module build definition

# Define the Geant4 Module.
geant4_add_module(G4hadronic_inclxx_interface
  PUBLIC_HEADERS
    G4INCLXXInterface.hh
    G4INCLXXInterfaceMessenger.hh
    G4INCLXXInterfaceStore.hh
    G4INCLXXVInterfaceTally.hh
  SOURCES
    G4INCLXXInterface.cc
    G4INCLXXInterfaceMessenger.cc
    G4INCLXXInterfaceStore.cc)

geant4_module_link_libraries(G4hadronic_inclxx_interface
  PUBLIC
    G4globman
    G4had_preequ_exciton
    G4hadronic_binary
    G4hadronic_deex_fission
    G4hadronic_deex_handler
    G4hadronic_deex_util
    G4hadronic_inclxx_physics
    G4hadronic_inclxx_utils
    G4hadronic_mgt
    G4hadronic_util
    G4intercoms
    G4partman
    G4track
  PRIVATE
    G4hadronic_abla
    G4hadronic_deex_evaporation
    G4hadronic_deex_management
    G4ions)
