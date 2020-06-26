#------------------------------------------------------------------------------
# Module : G4track
# Package: Geant4.src.G4track
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4track
  HEADERS
    G4FieldTrackUpdator.hh
    G4ForceCondition.hh
    G4GPILSelection.hh
    G4ParticleChange.hh
    G4ParticleChange.icc
    G4ParticleChangeForDecay.hh
    G4ParticleChangeForGamma.hh
    G4ParticleChangeForLoss.hh
    G4ParticleChangeForMSC.hh
    G4ParticleChangeForMSC.icc
    G4ParticleChangeForRadDecay.hh
    G4ParticleChangeForTransport.hh
    G4ParticleChangeForTransport.icc
    G4Step.hh
    G4Step.icc
    G4StepPoint.hh
    G4StepPoint.icc
    G4StepStatus.hh
    G4SteppingControl.hh
    G4Track.hh
    G4Track.icc
    G4TrackFastVector.hh
    G4TrackStatus.hh
    G4TrackVector.hh
    G4VParticleChange.hh
    G4VParticleChange.icc
    G4VelocityTable.hh
    G4VAuxiliaryTrackInformation.hh
    G4VUserTrackInformation.hh
    trkdefs.hh
  SOURCES
    G4FieldTrackUpdator.cc
    G4ParticleChange.cc
    G4ParticleChangeForDecay.cc
    G4ParticleChangeForGamma.cc
    G4ParticleChangeForLoss.cc
    G4ParticleChangeForMSC.cc
    G4ParticleChangeForTransport.cc
    G4Step.cc
    G4StepPoint.cc
    G4Track.cc
    G4VParticleChange.cc
    G4VelocityTable.cc
    G4VAuxiliaryTrackInformation.cc
    G4VUserTrackInformation.cc
  GRANULAR_DEPENDENCIES
    G4geometrymng
    G4globman
    G4intercoms
    G4magneticfield
    G4materials
    G4partman
  GLOBAL_DEPENDENCIES
    G4geometry
    G4global
    G4intercoms
    G4materials
    G4particles
  LINK_LIBRARIES
)

# List any source specific properties here
