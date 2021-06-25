# - G4UIGAG module build definition

# Define the Geant4 Module
geant4_add_module(G4UIGAG
  PUBLIC_HEADERS
    G4UIGAG.hh
    G4UIGainServer.hh
  SOURCES
    G4UIGAG.cc
    G4UIGainServer.cc)

geant4_module_link_libraries(G4UIGAG PUBLIC G4UIcommon G4intercoms G4globman)