#------------------------------------------------------------------------------
# Module : G4UIGAG
# Package: Geant4.src.G4interfaces.G4UIGAG
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4UIGAG
  HEADERS
    G4UIGAG.hh
    G4UIGainServer.hh
  SOURCES
    G4UIGAG.cc
    G4UIGainServer.cc
  GRANULAR_DEPENDENCIES
    G4UIcommon
    G4globman
    G4intercoms
  GLOBAL_DEPENDENCIES
    G4global
    G4intercoms
)

# List any source specific properties here
