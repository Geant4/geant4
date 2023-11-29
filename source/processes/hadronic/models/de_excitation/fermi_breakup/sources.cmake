# - G4hadronic_deex_fermi_breakup module build definition

# Define the Geant4 Module.
geant4_add_module(G4hadronic_deex_fermi_breakup
  PUBLIC_HEADERS
    G4FermiBreakUpVI.hh
    G4FermiChannels.hh
    G4FermiCoulombBarrier.hh
    G4FermiDecayProbability.hh
    G4FermiFragment.hh
    G4FermiFragmentsPoolVI.hh
    G4FermiPair.hh
    G4FermiPhaseSpaceDecay.hh
    G4VFermiBreakUp.hh
  SOURCES
    G4FermiBreakUpVI.cc
    G4FermiCoulombBarrier.cc
    G4FermiDecayProbability.cc
    G4FermiFragment.cc
    G4FermiFragmentsPoolVI.cc
    G4FermiPair.cc
    G4FermiPhaseSpaceDecay.cc)

geant4_module_link_libraries(G4hadronic_deex_fermi_breakup
  PUBLIC
    G4globman
    G4hadronic_deex_util
    G4hadronic_util
    G4hepgeometry
    G4heprandom
  PRIVATE
    G4hadronic_deex_management
    G4partman)
