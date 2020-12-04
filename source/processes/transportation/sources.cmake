#------------------------------------------------------------------------------
# sources.cmake
# Module : G4transportation
# Package: Geant4.src.G4processes.G4transportation
#
# Sources description for a library.
# Lists the sources and headers of the code explicitly.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
#
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
GEANT4_DEFINE_MODULE(NAME G4transportation
    HEADERS
        G4CoupledTransportation.hh
        G4CoupledTransportation.icc
        G4NeutronKiller.hh
        G4NeutronKillerMessenger.hh
        G4StepLimiter.hh
        G4TrackTerminator.hh
        G4Transportation.hh
        G4Transportation.icc
        G4UserSpecialCuts.hh
        G4VTrackTerminator.hh
	G4TransportationLogger.hh
	G4TransportationProcessType.hh
    SOURCES
        G4CoupledTransportation.cc
        G4NeutronKiller.cc
        G4NeutronKillerMessenger.cc
        G4StepLimiter.cc
        G4Transportation.cc
	G4TransportationLogger.cc
        G4UserSpecialCuts.cc
        G4VTrackTerminator.cc
    GRANULAR_DEPENDENCIES
        G4cuts
        G4emutils
        G4geombias
        G4geometrymng
        G4globman
        G4intercoms
        G4magneticfield
        G4materials
        G4navigation
        G4partman
        G4procman
        G4track
        G4volumes
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4intercoms
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

