#include "Tst50PhotonStandard.hh"

#include "G4ProcessManager.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

Tst50PhotonStandard::Tst50PhotonStandard(const G4String& name): G4VPhysicsConstructor(name)
{ }

Tst50PhotonStandard::~Tst50PhotonStandard()
{ }

void Tst50PhotonStandard::ConstructProcess()
{
  // Add standard processes for photons
  
  theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* manager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();

      if (particleName == "gamma") 
	{
          
	  manager->AddDiscreteProcess(new G4PhotoElectricEffect);
	  manager->AddDiscreteProcess(new G4ComptonScattering);
	  manager->AddDiscreteProcess(new G4GammaConversion);
          G4cout<<"Gamma Processes registered"<<G4endl;
	}   
    }
}
