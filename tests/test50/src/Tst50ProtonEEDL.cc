#include "Tst50ProtonEEDL.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4Proton.hh"
#include "G4hLowEnergyIonisation.hh"
#include "G4hLowEnergyLoss.hh"

Tst50ProtonEEDL::Tst50ProtonEEDL(const G4String& name): G4VPhysicsConstructor(name)
{ }

Tst50ProtonEEDL::~Tst50ProtonEEDL()
{ }

void Tst50ProtonEEDL::ConstructProcess()
{

theParticleIterator->reset();

  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* manager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      if (particleName == "proton") 
	{
	  // G4VProcess*  multipleScattering= new G4MultipleScattering(); 
   manager->AddProcess(new G4hLowEnergyIonisation(),-1,2,2);
   //  manager->AddProcess(multipleScattering,-1,1,1);  	
    G4hLowEnergyLoss::SetEnlossFluc(false);
	}}}
