#ifndef G4BCAction_h
#define G4BCAction_h 1
#include "G4CollisionInitialState.hh"
#include "G4KineticTrack.hh"
#include <vector>

class G4BCAction
{
  public:
  
  G4BCAction(){}
  virtual ~G4BCAction(){}
  virtual const std::vector<G4CollisionInitialState *> &
         GetCollisions(G4KineticTrack * aProjectile, 
	               std::vector<G4KineticTrack *> & someCandidates,
		       G4double aCurrentTime) = 0;

  virtual G4KineticTrackVector * 
         GetFinalState(G4KineticTrack * aProjectile, 
	               std::vector<G4KineticTrack *> & theTargets) = 0;
};
#endif
