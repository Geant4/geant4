#ifndef G4MesonAbsorption_hh
#define G4MesonAbsorption_hh

#include <vector>
#include "G4CollisionInitialState.hh"
#include "G4KineticTrack.hh"

class G4MesonAbsorption
{
  public:
  G4CollisionInitialState * GetCollision(G4KineticTrack * projectile, 
                                         vector<G4KineticTrack *> targets);
  vector<G4KineticTrack *> * Scatter(G4KineticTrack * projectile,  
                                     vector<G4KineticTrack *> & targets);
  private:
  
};

#endif
