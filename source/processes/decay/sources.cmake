# - G4decay module build definition

# Define the Geant4 Module.
geant4_add_module(G4decay
  PUBLIC_HEADERS
    G4Decay.hh
    G4DecayProcessType.hh
    G4DecayWithSpin.hh
    G4PionDecayMakeSpin.hh
    G4UnknownDecay.hh
    G4VExtDecayer.hh
  SOURCES
    G4Decay.cc
    G4DecayWithSpin.cc
    G4PionDecayMakeSpin.cc
    G4UnknownDecay.cc)

geant4_module_link_libraries(G4decay
  PUBLIC
    G4globman
    G4partman
    G4procman
    G4track
  PRIVATE
    G4hepgeometry
    G4heprandom
    G4magneticfield
    G4navigation)
