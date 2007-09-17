#ifndef G4GPRSELECTPARTICLES_HH
#define G4GPRSELECTPARTICLES_HH

#include "G4ParticleTable.hh"

class G4GPRSelectParticles {

public:

  template <typename Particle>
  void SelectParticle() 
  {
    fParticles.push_back(Particle::Definition());
  }

  void SelectParticle(G4ParticleDefinition* def) 
  {
    fParticles.push_back(def);
  }

  void SelectAllParticles() 
  { 
    G4ParticleTable::G4PTblDicIterator* iter = G4ParticleTable::GetParticleTable()->GetIterator();
    iter->reset();
    
    while ((*iter)()) {
      G4ParticleDefinition* particle = iter->value();
      fParticles.push_back(particle);
    }
  }
  /*
  void SelectParticleFunction(G4bool(*func)(G4ParticleDefinition*)) 
  {
    G4ParticleTable::G4PTblDicIterator* iter = G4ParticleTable::GetParticleTable()->GetIterator();
    iter->reset();
    
    while ((*iter)()) {      
      G4ParticleDefinition* particle = iter->value();
      if (func(particle)) fParticles.push_back(particle);
    }
  }
  */

  std::vector<G4ParticleDefinition*> GetList() const {return fParticles;}

private:

  std::vector<G4ParticleDefinition*> fParticles;

};

#endif
