#include "A01Biasing.hh"

#include "G4GPRBuilder.hh"
#include "G4GPRConverter.hh"
#include "G4Electron.hh"
#include "G4Transportation.hh"
#include "G4MultipleScattering.hh"
#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4GammaConversion.hh"
#include "G4ComptonScattering.hh"
#include "G4PhotoElectricEffect.hh"
#include "A01Triggers.hh"

/*

{
  G4cout<<"jane construct a01 biasing"<<G4endl;
  
  //  G4GPRConverter::LoadGPR(); // force generalised processing for all particles
  //G4GPRConverter::LoadGPR(G4Gamma::Definition()); // force generalised processing for all particles

   G4ProcessManager * pManager = 0;
   
   //Gamma
   pManager = G4Gamma::Gamma()->GetProcessManager();
   pManager->AddDiscreteProcess(new G4GammaConversion());
   pManager->AddDiscreteProcess(new G4ComptonScattering());
   pManager->AddDiscreteProcess(new G4PhotoElectricEffect());
   
   G4GPRConverter::LoadGPR(G4Gamma::Definition()); // force generalised processing for all particles
 
  // Override default A01 physics list for electrons. Just for fun, make the electrons use
  // the low energy physics list by default
//  G4String newPhysicsList("New Default");
//  CreateDefaultPhysicsList<G4Gamma>("New Default");
*/
using namespace G4GPRBuilder;
using namespace A01Triggers;

void A01Biasing::ConstructBiasing() 
{/*
  G4String caloListName("Calorimeter_Physics_List");

  CreatePhysicsListWithTrigger<G4Gamma, G4GPRTriggerTypes::Geometry::NewVolume>(caloListName, &CalorimeterTrigger);
  AddProcess<G4Gamma>(new G4Transportation,  -1, 0, 0, caloListName);

  CreatePhysicsListWithTrigger<G4Electron, G4GPRTriggerTypes::Geometry::NewVolume>(caloListName, &CalorimeterTrigger);
  AddProcess<G4Electron>(new G4Transportation,  -1, 0, 0, caloListName);

  CreatePhysicsListWithTrigger<G4Positron, G4GPRTriggerTypes::Geometry::NewVolume>(caloListName, &CalorimeterTrigger);
  AddProcess<G4Positron>(new G4Transportation,  -1, 0, 0, caloListName);
 */

  //  AddProcess<G4Gamma>(new G4GammaConversion, -1, -1, 1, caloListName);
  // AddProcess<G4Gamma>(new G4ComptonScattering, -1, -1, 2, caloListName);
  //AddProcess<G4Gamma>(new G4PhotoElectricEffect, -1, -1, 3, caloListName);
  /*
  AddProcess<G4Gamma>(new G4Transportation,  -1, 0, 0);
  AddProcess<G4Gamma>(new G4GammaConversion, -1, -1, 1);
  AddProcess<G4Gamma>(new G4ComptonScattering, -1, -1, 2);
  AddProcess<G4Gamma>(new G4PhotoElectricEffect, -1, -1, 3);
*/

  //  CreateDefaultPhysicsList<G4Electron>("New Default");

  /*
  AddProcess<G4Electron>(new G4Transportation,     -1, 0, 0);
  AddProcess<G4Electron>(new G4MultipleScattering, -1, 1, 1);
  AddProcess<G4Electron>(new G4eIonisation,        -1, 2, 2);
  AddProcess<G4Electron>(new G4eBremsstrahlung,    -1, 3, 3);
  */
  /*
  AddProcess<G4Electron>(new G4Transportation, -1, 0, 0);
  AddProcess<G4Electron>(new G4MultipleScattering,      -1, 1, 1);
  AddProcess<G4Electron>(new G4LowEnergyIonisation,     -1, 2, 2);
  AddProcess<G4Electron>(new G4LowEnergyBremsstrahlung, -1,-1, 3);
  */

  /*
 
  AddProcess<G4Positron>(new G4MultipleScattering, -1, 1, 1, newPhysicsList);
  AddProcess<G4Positron>(new G4eIonisation,        -1, 2, 2, newPhysicsList);
  AddProcess<G4Positron>(new G4eBremsstrahlung,    -1, 3, 3, newPhysicsList);
  AddProcess<G4Positron>(new G4eplusAnnihilation,   0,-1, 4, newPhysicsList);

  // Create yet another physics list for primary positrons
  CreatePhysicsListWithTrigger<G4Positron, G4GPRTriggerType::Tracking::StartTracking>("Primary Positrons", 
										      &A01GPRTriggers::PrimaryParticle);
  
  AddProcess<G4Positron>(new G4Transportation, ..);
  AddProcess<G4Positron>(new G4MultipleScattering);
  AddProcess<G4Positron>(new G4eplusAnnihilation);
    
  */

}
