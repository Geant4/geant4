#------------------------------------------------------------------------------
# Module : G4digits
# Package: Geant4.src.G4digits_hits.G4digits
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4digits
  HEADERS
    G4DCofThisEvent.hh
    G4TDigiCollection.hh
    G4VDigi.hh
    G4VDigiCollection.hh
  SOURCES
    G4DCofThisEvent.cc
    G4TDigiCollection.cc
    G4VDigi.cc
    G4VDigiCollection.cc
  GRANULAR_DEPENDENCIES
    G4globman
  GLOBAL_DEPENDENCIES
    G4global
)

# List any source specific properties here
