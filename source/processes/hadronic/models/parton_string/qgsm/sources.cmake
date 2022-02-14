# - G4hadronic_qgstring module build definition

# Define the Geant4 Module.
geant4_add_module(G4hadronic_qgstring
  PUBLIC_HEADERS
    G4BaryonSplitter.hh
    G4DiffractiveStringBuilder.hh
    G4GammaParticipants.hh
    G4MesonSplitter.hh
    G4PartonPair.hh
    G4QGSDiffractiveExcitation.hh
    G4QGSMSplitableHadron.hh
    G4QGSModel.hh
    G4QGSModel.icc
    G4QGSParticipants.hh
    G4QuarkExchange.hh
    G4Reggeons.hh
    G4SPBaryon.hh
    G4SPBaryonTable.hh
    G4SingleDiffractiveExcitation.hh
    G4SoftStringBuilder.hh
    G4SPPartonInfo.hh
  SOURCES
    G4BaryonSplitter.cc
    G4DiffractiveStringBuilder.cc
    G4GammaParticipants.cc
    G4MesonSplitter.cc
    G4PartonPair.cc
    G4QGSDiffractiveExcitation.cc
    G4QGSMSplitableHadron.cc
    G4QGSParticipants.cc
    G4QuarkExchange.cc
    G4Reggeons.cc
    G4SingleDiffractiveExcitation.cc
    G4SoftStringBuilder.cc
    G4SPBaryon.cc)

geant4_module_link_libraries(G4hadronic_qgstring
  PUBLIC
    G4baryons
    G4bosons
    G4globman
    G4had_string_man
    G4hadronic_util
    G4hepgeometry
    G4heprandom
    G4mesons
    G4partman)
