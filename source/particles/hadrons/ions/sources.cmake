# - G4ions module build definition

# Define the Geant4 Module.
geant4_add_module(G4ions
  PUBLIC_HEADERS
    G4Alpha.hh
    G4AntiAlpha.hh
    G4AntiDeuteron.hh
    G4AntiDoubleHyperDoubleNeutron.hh
    G4AntiDoubleHyperH4.hh
    G4AntiHe3.hh
    G4AntiHyperAlpha.hh
    G4AntiHyperHe5.hh
    G4AntiHyperH4.hh
    G4AntiHyperTriton.hh
    G4AntiTriton.hh
    G4Deuteron.hh
    G4DoubleHyperDoubleNeutron.hh
    G4DoubleHyperH4.hh
    G4GenericIon.hh
    G4He3.hh
    G4HyperAlpha.hh
    G4HyperHe5.hh
    G4HyperH4.hh
    G4HyperTriton.hh
    G4IonConstructor.hh
    G4Triton.hh
    G4GenericMuonicAtom.hh
  SOURCES
    G4Alpha.cc
    G4AntiAlpha.cc
    G4AntiDeuteron.cc
    G4AntiDoubleHyperDoubleNeutron.cc
    G4AntiDoubleHyperH4.cc
    G4AntiHe3.cc
    G4AntiHyperAlpha.cc
    G4AntiHyperHe5.cc
    G4AntiHyperH4.cc
    G4AntiHyperTriton.cc
    G4AntiTriton.cc
    G4Deuteron.cc
    G4DoubleHyperDoubleNeutron.cc
    G4DoubleHyperH4.cc
    G4GenericIon.cc
    G4He3.cc
    G4HyperAlpha.cc
    G4HyperHe5.cc
    G4HyperH4.cc
    G4HyperTriton.cc    
    G4IonConstructor.cc
    G4Triton.cc
    G4GenericMuonicAtom.cc)

geant4_module_link_libraries(G4ions PUBLIC G4globman G4partman)
