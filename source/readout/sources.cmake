# - G4readout module build definition

# Define the Geant4 Module.
geant4_add_module(G4readout
  PUBLIC_HEADERS
    G4DCtable.hh
    G4DMmessenger.hh
    G4DigiManager.hh
    G4VDigitizerModule.hh
  SOURCES
    G4DCtable.cc
    G4DMmessenger.cc
    G4DigiManager.cc
    G4VDigitizerModule.cc)

geant4_module_link_libraries(G4readout
  PUBLIC
    G4globman
    G4intercoms
  PRIVATE
    G4detector
    G4digits
    G4event
    G4hits
    G4run)
