
#include "Particles.hh"
#include "G4LeptonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"


Particles::Particles(const G4String& name) :  
   G4VPhysicsConstructor(name) { 

}


Particles::~Particles() {

}


void Particles::ConstructParticle() { 

  G4LeptonConstructor lepton;
  lepton.ConstructParticle();
 
  G4BosonConstructor boson;
  boson.ConstructParticle();

  G4MesonConstructor meson;
  meson.ConstructParticle();

  G4BaryonConstructor baryon;
  baryon.ConstructParticle();

  G4ShortLivedConstructor shortLived;
  shortLived.ConstructParticle();
}
