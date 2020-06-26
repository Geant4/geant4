#------------------------------------------------------------------------------
# Module : G4readout
# Package: Geant4.src.G4readout
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4readout
  HEADERS
    G4DCtable.hh
    G4DMmessenger.hh
    G4DigiManager.hh
  G4VDigitizerModule.hh
  SOURCES
    G4DCtable.cc
    G4DMmessenger.cc
    G4DigiManager.cc
    G4VDigitizerModule.cc
  GRANULAR_DEPENDENCIES
    G4detector
    G4digits
    G4emutils
    G4event
    G4geometrymng
    G4globman
    G4hits
    G4intercoms
    G4materials
    G4navigation
    G4partman
    G4procman
    G4run
    G4track
    G4tracking
    G4volumes
  GLOBAL_DEPENDENCIES
    G4digits_hits
    G4event
    G4geometry
    G4global
    G4intercoms
    G4materials
    G4particles
    G4processes
    G4run
    G4track
    G4tracking
  LINK_LIBRARIES
)

# List any source specific properties here

