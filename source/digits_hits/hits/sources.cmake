#------------------------------------------------------------------------------
# Module : G4hits
# Package: Geant4.src.G4digits_hits.G4hits
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4hits
  HEADERS
    G4HCofThisEvent.hh
    G4THitsCollection.hh
    G4THitsMap.hh
    G4THitsVector.hh
    G4VHit.hh
    G4VHitsCollection.hh
  SOURCES
    G4HCofThisEvent.cc
    G4THitsCollection.cc
    G4VHit.cc
    G4VHitsCollection.cc
  GRANULAR_DEPENDENCIES
    G4globman
  GLOBAL_DEPENDENCIES
    G4global
)

# List any source specific properties here
