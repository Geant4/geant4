#include "PhysListEmLivermore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyRayleigh.hh"

#include "G4MultipleScattering.hh"

#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmLivermore::PhysListEmLivermore(const G4String& name)
  :  G4VPhysicsConstructor(name)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmLivermore::~PhysListEmLivermore()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmLivermore::ConstructProcess()
{
  // Add LowEn or standard EM Processes

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
      // gamma         
      pmanager->AddDiscreteProcess(new G4LowEnergyPhotoElectric);
      pmanager->AddDiscreteProcess(new G4LowEnergyCompton);
      pmanager->AddDiscreteProcess(new G4LowEnergyGammaConversion);
      pmanager->AddDiscreteProcess(new G4LowEnergyRayleigh);
      
    } else if (particleName == "e-") {
      //electron
      pmanager->AddProcess(new G4MultipleScattering,      -1, 1, 1);
      pmanager->AddProcess(new G4LowEnergyIonisation,     -1, 2, 2);
      pmanager->AddProcess(new G4LowEnergyBremsstrahlung, -1,-1, 3);
	    
    } else if (particleName == "e+") {
      //positron
      pmanager->AddProcess(new G4MultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4eIonisation,        -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1, 3, 3);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1, 4);
      
    } 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

