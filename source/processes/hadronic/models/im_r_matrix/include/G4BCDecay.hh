#ifndef G4BCDecay_h
#define G4BCDecay_h 1
#include "G4BCAction.hh"
#include "G4CollisionInitialState.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include <vector>

class G4BCDecay : public G4BCAction
{
  public:
  
  G4BCDecay(){}
  virtual ~G4BCDecay(){}
  virtual const std::vector<G4CollisionInitialState *> &
         GetCollisions(G4KineticTrack * aProjectile, 
	               std::vector<G4KineticTrack *> & ,
		       G4double theCurrentTime)
  {
    theColl.clear();
    if(aProjectile->GetDefinition()->IsShortLived())
    {
      G4double aTime = theCurrentTime+aProjectile->SampleResidualLifetime();
      G4KineticTrackVector noTarget;
      G4CollisionInitialState * aDecay = 
            new G4CollisionInitialState(aTime, aProjectile, noTarget, this);
      theColl.push_back(aDecay);
    }
    return theColl;
  }

  virtual G4KineticTrackVector * GetFinalState(G4KineticTrack * aProjectile, 
	                           std::vector<G4KineticTrack *> & )
  {
    return aProjectile->Decay();
  }
  
  private:
  
  std::vector<G4CollisionInitialState *> theColl;
};

#endif
