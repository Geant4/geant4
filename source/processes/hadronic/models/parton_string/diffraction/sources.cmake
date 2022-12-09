# - G4had_string_diff module build definition

# Define the Geant4 Module.
geant4_add_module(G4had_string_diff
  PUBLIC_HEADERS
    G4DiffractiveExcitation.hh
    G4DiffractiveSplitableHadron.hh
    G4ElasticHNScattering.hh
    G4FTFAnnihilation.hh
    G4FTFModel.hh
    G4FTFParameters.hh
    G4FTFParticipants.hh
    G4FTFTunings.hh
    G4FTFTuningsMessenger.hh
  SOURCES
    G4DiffractiveExcitation.cc
    G4DiffractiveSplitableHadron.cc
    G4ElasticHNScattering.cc
    G4FTFAnnihilation.cc
    G4FTFModel.cc
    G4FTFParameters.cc
    G4FTFParticipants.cc
    G4FTFTunings.cc
    G4FTFTuningsMessenger.cc)

geant4_module_link_libraries(G4had_string_diff
  PUBLIC
    G4baryons
    G4globman
    G4had_string_man
    G4hadronic_util
    G4hepgeometry
    G4intercoms
  PRIVATE
    G4had_string_frag
    G4hadronic_xsect
    G4heprandom
    G4mesons
    G4partman)
