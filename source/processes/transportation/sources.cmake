# - G4transportation module build definition

# Define the Geant4 Module.
geant4_add_module(G4transportation
  PUBLIC_HEADERS
    G4CoupledTransportation.hh
    G4NeutronKiller.hh
    G4NeutronKillerMessenger.hh
    G4StepLimiter.hh
    G4TrackTerminator.hh
    G4Transportation.hh
    G4Transportation.icc
    G4VTrackTerminator.hh
    G4TransportationLogger.hh
    G4TransportationParameters.hh
    G4TransportationProcessType.hh
  SOURCES
    G4CoupledTransportation.cc
    G4NeutronKiller.cc
    G4NeutronKillerMessenger.cc
    G4StepLimiter.cc
    G4Transportation.cc
    G4TransportationParameters.cc    
    G4TransportationLogger.cc
    G4VTrackTerminator.cc)

geant4_module_link_libraries(G4transportation
  PUBLIC
    G4globman
    G4intercoms
    G4partman
    G4procman
    G4track
  PRIVATE
    G4cuts
    G4magneticfield
    G4navigation)
