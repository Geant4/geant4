// $Id: ExE03PhysicsList.cc,v 1.2 1999/04/17 05:27:14 kurasige Exp $

#include "globals.hh"
#include "ExE03PhysicsList.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BarionConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include <iomanip.h>                


ExE03PhysicsList::ExE03PhysicsList():  G4VUserPhysicsList()
{
  SetVerboseLevel(1);
}

ExE03PhysicsList::~ExE03PhysicsList()
{
}

void ExE03PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBarions();
  ConstructIons();
}

void ExE03PhysicsList::ConstructBosons()
{
  // Construct all bosons
  G4BosonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void ExE03PhysicsList::ConstructLeptons()
{
  // Construct all leptons
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void ExE03PhysicsList::ConstructMesons()
{
  //  Construct all mesons
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void ExE03PhysicsList::ConstructBarions()
{
  //  Construct all barions
  G4BarionConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void ExE03PhysicsList::ConstructIons()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}


void ExE03PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructLeptHad();
  ConstructHad();
  ConstructGeneral();
}

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
void ExE03PhysicsList::ConstructEM()
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
  
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
    //muon  
     // Construct processes for muon+
     pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);
     pmanager->AddProcess(new G4MuIonisation(),-1,2,2);
     pmanager->AddProcess(new G4MuBremsstrahlung(),-1,-1,3);
     pmanager->AddProcess(new G4MuPairProduction(),-1,-1,4);       
     
    } else { 
      if ((particle->GetPDGCharge() != 0.0) && 
          (particle->GetParticleName() != "chargedgeantino")) {
     // all others charged particles except geantino
       pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);
       pmanager->AddProcess(new G4hIonisation(),-1,2,2);       
     }
    }
  }
}

void ExE03PhysicsList::ConstructHad()
{;}

void ExE03PhysicsList::ConstructLeptHad()
{;}

#include "G4Decay.hh"
void ExE03PhysicsList::ConstructGeneral()
{
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager->AddProcess(theDecayProcess, INT_MAX, -1, INT_MAX); 
   }
  }
}

void ExE03PhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "ExE03PhysicsList::SetCuts:";
    G4cout << "CutLength : " << defaultCutValue/mm << " (mm)" << endl;
  }  

  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}
 
