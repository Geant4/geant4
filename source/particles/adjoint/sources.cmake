# - G4partadj module build definition

# Define the Geant4 Module.
geant4_add_module(G4partadj
  PUBLIC_HEADERS
    G4AdjointAlpha.hh
    G4AdjointDeuteron.hh
    G4AdjointElectron.hh
    G4AdjointElectronFI.hh
    G4AdjointGamma.hh
    G4AdjointGenericIon.hh
    G4AdjointHe3.hh
    G4AdjointIons.hh
    G4AdjointPositron.hh
    G4AdjointProton.hh
    G4AdjointTriton.hh
  SOURCES
    G4AdjointAlpha.cc
    G4AdjointDeuteron.cc
    G4AdjointElectron.cc
    G4AdjointElectronFI.cc
    G4AdjointGamma.cc
    G4AdjointGenericIon.cc
    G4AdjointHe3.cc
    G4AdjointIons.cc
    G4AdjointPositron.cc
    G4AdjointProton.cc
    G4AdjointTriton.cc)

geant4_module_link_libraries(G4partadj PUBLIC G4partman G4globman)
