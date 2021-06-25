# - G4shortlived module build definition

# Define the Geant4 Module.
geant4_add_module(G4shortlived
  PUBLIC_HEADERS
    G4DiQuarks.hh
    G4ExcitedBaryonConstructor.hh
    G4ExcitedBaryons.hh
    G4ExcitedDeltaConstructor.hh
    G4ExcitedLambdaConstructor.hh
    G4ExcitedMesonConstructor.hh
    G4ExcitedMesons.hh
    G4ExcitedNucleonConstructor.hh
    G4ExcitedSigmaConstructor.hh
    G4ExcitedXiConstructor.hh
    G4Gluons.hh
    G4Quarks.hh
    G4ShortLivedConstructor.hh
    G4VShortLivedParticle.hh
  SOURCES
    G4DiQuarks.cc
    G4ExcitedBaryonConstructor.cc
    G4ExcitedBaryons.cc
    G4ExcitedDeltaConstructor.cc
    G4ExcitedLambdaConstructor.cc
    G4ExcitedMesonConstructor.cc
    G4ExcitedMesons.cc
    G4ExcitedNucleonConstructor.cc
    G4ExcitedSigmaConstructor.cc
    G4ExcitedXiConstructor.cc
    G4Gluons.cc
    G4Quarks.cc
    G4ShortLivedConstructor.cc
    G4VShortLivedParticle.cc)

geant4_module_link_libraries(G4shortlived PUBLIC G4globman G4partman)
