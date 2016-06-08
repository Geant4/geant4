// $Id: ExE01PhysicsList.cc,v 1.2 1999/04/17 04:06:46 kurasige Exp $

#include "globals.hh"
#include "ExE01PhysicsList.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include <iomanip.h>                


ExE01PhysicsList::ExE01PhysicsList():  G4VUserPhysicsList()
{
  SetVerboseLevel(1);
}

ExE01PhysicsList::~ExE01PhysicsList()
{
}

void ExE01PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();

}

void ExE01PhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}

void ExE01PhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
}

void ExE01PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
}

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

void ExE01PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
    // gamma
      // Construct processes for gamma
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());      
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

    } else if (particleName == "e-") {
    //electron
      // Construct processes for electron
      pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);
      pmanager->AddProcess(new G4eIonisation(),-1,2,2);
      pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);
  
    } else if (particleName == "e+") {
    //positron
      // Construct processes for positron
     pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);
     
     pmanager->AddProcess(new G4eIonisation(),-1,2,2);
     pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);      
     pmanager->AddProcess(new G4eplusAnnihilation(),0,-1,4);
  
    }
  }
}

void ExE01PhysicsList::SetCuts()
{
  G4double cut = defaultCutValue;
  if (verboseLevel >0){
    G4cout << "ExE01PhysicsList::SetCuts:";
    G4cout << "CutLength : " << cut/mm << " (mm)" << endl;
  }  

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  SetCutValue(cut, "gamma");
  SetCutValue(cut, "e-");
  SetCutValue(cut, "e+");
 
  SetCutValueForOthers(cut);

  if (verboseLevel>1) {
    DumpCutValuesTable();
  }
}
