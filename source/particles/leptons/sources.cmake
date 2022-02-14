# - G4leptons module build definition

# Define the Geant4 Module.
geant4_add_module(G4leptons
  PUBLIC_HEADERS
    G4AntiNeutrinoE.hh
    G4AntiNeutrinoMu.hh
    G4AntiNeutrinoTau.hh
    G4Electron.hh
    G4LeptonConstructor.hh
    G4MuonMinus.hh
    G4MuonPlus.hh
    G4NeutrinoE.hh
    G4NeutrinoMu.hh
    G4NeutrinoTau.hh
    G4Positron.hh
    G4TauMinus.hh
    G4TauPlus.hh
  SOURCES
    G4AntiNeutrinoE.cc
    G4AntiNeutrinoMu.cc
    G4AntiNeutrinoTau.cc
    G4Electron.cc
    G4LeptonConstructor.cc
    G4MuonMinus.cc
    G4MuonPlus.cc
    G4NeutrinoE.cc
    G4NeutrinoMu.cc
    G4NeutrinoTau.cc
    G4Positron.cc
    G4TauMinus.cc
    G4TauPlus.cc)

geant4_module_link_libraries(G4leptons PUBLIC G4globman G4partman)
