#include "G4StringChipsInterface.hh"
#include "global.hh"
#include "G4Pair.hh"

G4VParticleChange* G4StringChipsInterface::
ApplyYourself(const G4Track& aTrack, G4Nucleus& theNucleus)
{
  G4Exception("G4StringChipsInterface: No stand-alone mode for all particles, please use specialised model classes");
  return NULL;
}
G4ReactionProductVector* G4StringChipsInterface::
Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus)
{
  G4int aTrack;
  if(theSecondaries->length() == 1) 
    G4Exception("G4StringChipsInterface: Only one particle from String models!");
  G4Pair<G4double, G4double> theImpact = theNucleus->RefetchImpactXandY();
  G4double impactX = theImpact.first;
  G4double impactY = theImpact.second;
  G4double inpactPar2 = impactX*impactX + impactY*impactY;
  
  G4double radius2 = theNucleus->GetNuclearRadius(5*percent);
  G4double pathlength = 2.*sqrt(radius2 - inpactPar2);
  G4double energyLostInFragmentation = 1*GeV*pathlength/fermi;
  
  // now select all particles in range
  
  // now call chips with this info in place
} 
