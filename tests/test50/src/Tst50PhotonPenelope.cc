
#include "Tst50PhotonPenelope.hh"

#include "G4ProcessManager.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4PenelopeCompton.hh"
#include "G4PenelopeGammaConversion.hh"
#include "G4PenelopePhotoElectric.hh"
#include "G4PenelopeRayleigh.hh"

Tst50PhotonPenelope::Tst50PhotonPenelope(const G4String& name): G4VPhysicsConstructor(name)
{ }

Tst50PhotonPenelope::~Tst50PhotonPenelope()
{ }

void Tst50PhotonPenelope::ConstructProcess()
{
  
  
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* manager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      if (particleName == "gamma") 
	{
	  manager->AddDiscreteProcess(new G4PenelopePhotoElectric);
	  manager->AddDiscreteProcess(new G4PenelopeCompton);
	  manager->AddDiscreteProcess(new G4PenelopeGammaConversion);
	  manager->AddDiscreteProcess(new G4PenelopeRayleigh);
	}   
    }
}


