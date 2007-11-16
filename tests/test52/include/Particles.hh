#ifndef PARTICLES_HH
#define PARTICLES_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"


class Particles : public G4VPhysicsConstructor {

 public: 
   Particles(const G4String& name = "Particles");
   virtual ~Particles();
 
 protected: 
   virtual void ConstructParticle();
   virtual void ConstructProcess() {}; 
};

#endif // PARTICLES_HH








