#ifndef G4ParticleChangeForNothing_hh
#define G4ParticleChangeForNothing_hh 1

#include "G4VParticleChange.hh"

class G4ParticleChangeForNothing : public G4VParticleChange {
public:
  G4ParticleChangeForNothing() : G4VParticleChange() {}
  ~G4ParticleChangeForNothing() {}

public:
  // -- from base class G4VParticleChange:
  virtual void Initialize(const G4Track &track)
  {
    theStatusChange = track.GetTrackStatus();
    theNumberOfSecondaries = 0;
  }
  virtual G4Step* UpdateStepForAtRest   (G4Step* step) {return step;}
  virtual G4Step* UpdateStepForAlongStep(G4Step* step) {return step;}
  virtual G4Step* UpdateStepForPostStep (G4Step* step) {return step;}
};

#endif
