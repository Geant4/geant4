# - G4hits module build definition

# Define the Geant4 Module.
geant4_add_module(G4hits
  PUBLIC_HEADERS
    G4HCofThisEvent.hh
    G4THitsCollection.hh
    G4THitsMap.hh
    G4THitsVector.hh
    G4VHit.hh
    G4VHitsCollection.hh
  SOURCES
    G4HCofThisEvent.cc
    G4THitsCollection.cc
    G4VHitsCollection.cc)

geant4_module_link_libraries(G4hits PUBLIC G4globman)
