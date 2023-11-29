# - G4physlist_util module build definition

# Define the Geant4 Module.
geant4_add_module(G4physlist_util
  PUBLIC_HEADERS
    CompileTimeConstraints.hh
    G4HadParticles.hh
    G4HadProcesses.hh
    G4PhysListUtil.hh
    G4WarnPLStatus.hh
  SOURCES
    G4HadParticles.cc
    G4HadProcesses.cc
    G4PhysListUtil.cc
    G4WarnPLStatus.cc)

geant4_module_link_libraries(G4physlist_util
  PUBLIC
    G4globman
    G4hadronic_mgt
    G4hadronic_xsect
    G4hadronic_proc
    G4partman
    G4procman
  PRIVATE
    G4baryons
    G4emutils
    G4hadronic_deex_management
    G4hadronic_deex_photon_evaporation
    G4hadronic_util
    G4had_string_diff
    G4materials
    G4run)
