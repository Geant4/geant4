#ifndef G4COLLIDER_HH
#define G4COLLIDER_HH

//#ifndef G4INUCL_PARTICLE_HH
#include "G4InuclParticle.hh"
//#endif

#include "G4CollisionOutput.hh"

class G4Collider {

public:

  G4Collider() {};

  virtual G4CollisionOutput collide(G4InuclParticle* bullet,
				    G4InuclParticle* target) = 0;

};        

#endif // G4COLLIDER_HH 
