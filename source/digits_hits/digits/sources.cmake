# - G4digits module build definition

# Define the Geant4 Module.
geant4_add_module(G4digits
  PUBLIC_HEADERS
    G4DCofThisEvent.hh
    G4TDigiCollection.hh
    G4VDigi.hh
    G4VDigiCollection.hh
  SOURCES
    G4DCofThisEvent.cc
    G4TDigiCollection.cc
    G4VDigi.cc
    G4VDigiCollection.cc)

geant4_module_link_libraries(G4digits PUBLIC G4globman)
