# - G4had_string_man module build definition

# Define the Geant4 Module.
geant4_add_module(G4had_string_man
  PUBLIC_HEADERS
    G4InteractionContent.hh
    G4StringModel.hh
    G4VParticipants.hh
    G4VPartonStringModel.hh
    G4VSplitableHadron.hh
    G4VStringFragmentation.hh
  SOURCES
    G4InteractionContent.cc
    G4StringModel.cc
    G4VParticipants.cc
    G4VPartonStringModel.cc
    G4VSplitableHadron.cc
    G4VStringFragmentation.cc)

geant4_module_link_libraries(G4had_string_man
  PUBLIC
    G4globman
    G4hadronic_mgt
    G4hadronic_util
    G4hepgeometry
    G4partman
  PRIVATE
    G4shortlived)
