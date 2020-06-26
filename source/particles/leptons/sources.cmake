#------------------------------------------------------------------------------
# Module : G4leptons
# Package: Geant4.src.G4particles.G4leptons
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4leptons
  HEADERS
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
    G4TauPlus.cc
  GRANULAR_DEPENDENCIES
    G4globman
    G4materials
    G4partman
  GLOBAL_DEPENDENCIES
    G4global
    G4materials
)

# List any source specific properties here
