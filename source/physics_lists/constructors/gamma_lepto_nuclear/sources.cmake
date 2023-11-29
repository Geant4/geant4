# - G4phys_ctor_glnuclear module build definition

# Define the Geant4 Module.
geant4_add_module(G4phys_ctor_glnuclear
  PUBLIC_HEADERS
    G4BertiniElectroNuclearBuilder.hh
    G4EmExtraPhysics.hh
    G4EmMessenger.hh
    G4LENDBertiniGammaElectroNuclearBuilder.hh
  SOURCES
    G4BertiniElectroNuclearBuilder.cc
    G4EmExtraPhysics.cc
    G4EmMessenger.cc
    G4LENDBertiniGammaElectroNuclearBuilder.cc)

geant4_module_link_libraries(G4phys_ctor_glnuclear
  PUBLIC
    G4globman
    G4had_lept_nuclear
    G4had_string_frag
    G4had_theo_max
    G4hadronic_bert_cascade
    G4hadronic_binary
    G4hadronic_proc
    G4hadronic_qgstring
    G4intercoms
    G4run
  PRIVATE
    G4baryons
    G4bosons
    G4emhighenergy
    G4emutils
    G4had_gamma_nuclear
    G4had_lend
    G4had_preequ_exciton
    G4hadronic_coherent_elastic
    G4hadronic_util
    G4hadronic_xsect
    G4ions
    G4leptons
    G4mesons
    G4muons
    G4partman
    G4phys_builders
    G4phys_ctor_em
    G4phys_ctor_factory
    G4procman
    G4xrays)
