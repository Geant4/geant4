#include "Tst50ProtonEEDLziegler.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4MultipleScattering.hh"
#include "G4Proton.hh"
#include "G4hLowEnergyIonisation.hh"
#include "G4hLowEnergyLoss.hh"
#include "G4hZiegler1985p.hh"

Tst50ProtonEEDLziegler::Tst50ProtonEEDLziegler(const G4String& name): G4VPhysicsConstructor(name)
{ }

Tst50ProtonEEDLziegler::~Tst50ProtonEEDLziegler()
{ }

void Tst50ProtonEEDLziegler::ConstructProcess()
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
	  G4hLowEnergyIonisation* ion= new G4hLowEnergyIonisation();
   manager->AddProcess(ion,-1,2,2);

   ion ->SetElectronicStoppingPowerModel(particle, "G4hZiegler1985p");
  
 //  manager->AddProcess(multipleScattering,-1,1,1);  	
    G4hLowEnergyLoss::SetEnlossFluc(false);
	}}}
