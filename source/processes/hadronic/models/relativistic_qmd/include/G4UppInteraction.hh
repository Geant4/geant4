
#ifndef G4UPPINTERACTION_H
#define G4UPPINTERACTION_H


#include "G4UppTrackVector.hh"
#include "G4UppTrackChange.hh"
#include "g4std/vector"


class G4UppInteraction
{
public:

  G4UppInteraction() : InteractionTime(-1) {}
  G4UppInteraction(const G4double Time, const G4int i, const G4int j);

  G4UppTrackChange Perform();

  G4double getInteractionTime() 
    { return InteractionTime; }
  G4UppTrackVector* getOutgoingParticles() 
    { return &OutgoingParticles; }

private:

  G4std::vector<G4int> IncomingParticleIndex;
  G4UppTrackVector OutgoingParticles;
  G4double InteractionTime;

};


#endif // G4UPPINTERACTION_H



