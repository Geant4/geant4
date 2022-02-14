# - G4emdna-moltypes model build definition

# Define the Geant4 Module.
geant4_add_module(G4emdna-moltypes
  PUBLIC_HEADERS
    G4Electron_aq.hh
    G4FakeMolecule.hh
    G4H2.hh
    G4H2O2.hh
    G4H2O.hh
    G4H3O.hh
    G4HO2.hh
    G4Hydrogen.hh
    G4O2.hh
    G4O3.hh
    G4OH.hh
    G4Oxygen.hh
    G4DNAMolecule.hh
  SOURCES
    G4Electron_aq.cc
    G4FakeMolecule.cc
    G4H2.cc
    G4H2O2.cc
    G4H2O.cc
    G4H3O.cc
    G4HO2.cc
    G4Hydrogen.cc
    G4O2.cc
    G4O3.cc
    G4OH.cc
    G4Oxygen.cc
    G4DNAMolecule.cc)

geant4_module_link_libraries(G4emdna-moltypes
  PUBLIC
    G4emdna-molman
    G4globman
    G4partman)
