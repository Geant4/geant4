#------------------------------------------------------------------------------
# sources.cmake
# Module : G4biasing_gen
# Package: Geant4.src.G4processes.G4biasing_gen
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
GEANT4_DEFINE_MODULE(NAME G4biasing_gen
    HEADERS
        G4BiasingHelper.hh
        G4BiasingProcessInterface.hh
        G4BiasingProcessSharedData.hh
        G4BOptnChangeCrossSection.hh
        G4BOptnCloning.hh
        G4BOptnForceCommonTruncatedExp.hh
        G4BOptnForceFreeFlight.hh
	G4BOptnLeadingParticle.hh
        G4BOptrForceCollision.hh
        G4BOptrForceCollisionTrackData.hh
        G4ILawCommonTruncatedExp.hh
        G4ILawForceFreeFlight.hh
        G4ILawTruncatedExp.hh
        G4InteractionLawPhysical.hh
	G4ParallelGeometriesLimiterProcess.hh
        G4ParticleChangeForNothing.hh
        G4ParticleChangeForOccurenceBiasing.hh
    SOURCES
        G4BiasingHelper.cc
        G4BiasingProcessInterface.cc
	G4BiasingProcessSharedData.cc
        G4BOptnChangeCrossSection.cc
        G4BOptnCloning.cc
        G4BOptnForceCommonTruncatedExp.cc
        G4BOptnForceFreeFlight.cc
	G4BOptnLeadingParticle.cc
        G4BOptrForceCollision.cc
	G4BOptrForceCollisionTrackData.cc
        G4ILawCommonTruncatedExp.cc
        G4ILawForceFreeFlight.cc
        G4ILawTruncatedExp.cc
        G4InteractionLawPhysical.cc
	G4ParallelGeometriesLimiterProcess.cc
        G4ParticleChangeForOccurenceBiasing.cc
    GRANULAR_DEPENDENCIES
        G4cuts
        G4detector
        G4emutils
        G4geombias
        G4geometrymng
        G4globman
        G4hits
        G4intercoms
        G4magneticfield
        G4materials
        G4navigation
        G4partman
        G4procman
        G4track
        G4transportation
        G4volumes
	G4biasing_mgt
    GLOBAL_DEPENDENCIES
        G4digits_hits
        G4geometry
        G4global
        G4intercoms
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

