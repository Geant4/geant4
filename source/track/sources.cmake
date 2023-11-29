# - G4track module build definition

# Define the Geant4 Module.
geant4_add_module(G4track
  PUBLIC_HEADERS
    G4FieldTrackUpdator.hh
    G4ForceCondition.hh
    G4GPILSelection.hh
    G4ParticleChange.hh
    G4ParticleChange.icc
    G4ParticleChangeForDecay.hh
    G4ParticleChangeForGamma.hh
    G4ParticleChangeForLoss.hh
    G4ParticleChangeForMSC.hh
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
    G4VUserTrackInformation.cc)

geant4_module_compile_definitions(G4track PRIVATE G4TRACK_ALLOC_EXPORT)

geant4_module_link_libraries(G4track
  PUBLIC
    G4geometrymng
    G4globman
    G4materials
    G4partman
  PRIVATE
    G4magneticfield)
