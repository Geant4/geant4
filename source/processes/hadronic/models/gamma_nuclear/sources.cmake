# - G4had_gamma_nuclear module build definition

# Define the Geant4 Module.
geant4_add_module(G4had_gamma_nuclear
  PUBLIC_HEADERS G4LENDorBERTModel.hh
  SOURCES G4LENDorBERTModel.cc)

geant4_module_link_libraries(G4had_gamma_nuclear
  PUBLIC
    G4had_lend
  PRIVATE
    G4globman
    G4hadronic_bert_cascade
    G4partman)
