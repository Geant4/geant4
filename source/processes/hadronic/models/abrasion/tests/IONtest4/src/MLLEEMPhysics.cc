////////////////////////////////////////////////////////////////////////////////
//
#include "MLLEEMPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
////////////////////////////////////////////////////////////////////////////////
//
MLLEEMPhysics::MLLEEMPhysics (const G4String& name)
  :G4VPhysicsConstructor(name)
{}
////////////////////////////////////////////////////////////////////////////////
//
MLLEEMPhysics::~MLLEEMPhysics ()
{}
////////////////////////////////////////////////////////////////////////////////
//
#include "G4ParticleDefinition.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleTable.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"

void MLLEEMPhysics::ConstructParticle ()
{
  // gamma
  G4Gamma::GammaDefinition();
 
  // electron
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
}
////////////////////////////////////////////////////////////////////////////////
//
#include "G4ProcessManager.hh"

// Electromagnetic Processes ////////////////////////////////////////////////
// all charged particles
#include "G4MultipleScattering.hh"

// gamma
#include "G4LowEnergyRayleigh.hh" 
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyCompton.hh"  
#include "G4LowEnergyGammaConversion.hh" 


// e-
#include "G4LowEnergyIonisation.hh" 
#include "G4LowEnergyBremsstrahlung.hh" 

// e+
#include "G4eIonisation.hh" 
#include "G4eBremsstrahlung.hh" 
#include "G4eplusAnnihilation.hh"

void MLLEEMPhysics::ConstructProcess ()
{
  theParticleIterator->reset();
  //processes
  G4LowEnergyPhotoElectric* lowePhot = new G4LowEnergyPhotoElectric();
  G4LowEnergyIonisation* loweIon  = new G4LowEnergyIonisation();
  G4LowEnergyBremsstrahlung* loweBrem = new G4LowEnergyBremsstrahlung();
  
  // note LowEIon uses proton as basis for its data-base, therefore
  // cannot specify different LowEnergyIonisation models for different
  // particles, but can change model globally for Ion, Alpha and Proton.
  G4MultipleScattering* aMultipleScattering = new G4MultipleScattering();
  
  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    G4String particleType = particle->GetParticleType();
    //    G4double charge = particle->GetPDGCharge();
    if (particleName == "gamma") {
      //gamma
      pmanager->AddDiscreteProcess(new G4LowEnergyRayleigh());
      pmanager->AddDiscreteProcess(lowePhot);
      pmanager->AddDiscreteProcess(new G4LowEnergyCompton());
      pmanager->AddDiscreteProcess(new G4LowEnergyGammaConversion());

    } else if (particleName == "e-") {
      //electron
      // process ordering: AddProcess(name, at rest, along step, post step)
      // -1 = not implemented, then ordering............
      pmanager->AddProcess(aMultipleScattering,    -1, 1,1);
      pmanager->AddProcess(loweIon,                       -1, 2,2);
      pmanager->AddProcess(loweBrem,                      -1,-1,3);

    } else if (particleName == "e+") {
    //positron
      pmanager->AddProcess(aMultipleScattering, -1, 1,1);
      pmanager->AddProcess(new G4eIonisation(),        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung(),    -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation(),   0,-1,4);      
      
    }

    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(250*eV,1e5*MeV);
    //
    //fluorescence apply specific cut for flourescence from photons, electrons
    //and bremsstrahlung photons:
    G4double cut = 250*eV;
    lowePhot->SetCutForLowEnSecPhotons(cut);
    loweIon->SetCutForLowEnSecPhotons(cut);
    loweBrem->SetCutForLowEnSecPhotons(cut);
  }
}
////////////////////////////////////////////////////////////////////////////////
