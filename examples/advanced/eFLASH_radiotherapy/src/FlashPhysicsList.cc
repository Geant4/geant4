//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file FlashPhysicsList.cc
/// \brief Implementation of the FlashPhysicsList class

#include "FlashPhysicsList.hh"
#include "FlashPhysicsListMessenger.hh"
#include "G4DecayPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4PhysListFactory.hh"
#include "G4ProductionCuts.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4SystemOfUnits.hh"

FlashPhysicsList::FlashPhysicsList() : G4VModularPhysicsList() {



    pMessenger = new FlashPhysicsListMessenger(this);
      SetVerboseLevel(2);
  
// EM physics
  emPhysicsList = new G4EmStandardPhysics_option4(1);
  emName = G4String("emstandard_opt4");

   // Decay physics and all particles
  decPhysicsList = new G4DecayPhysics();

     // Radioactive physics and all particles
  radPhysicsList=new G4RadioactiveDecayPhysics();


}

FlashPhysicsList::~FlashPhysicsList() {

delete pMessenger;
  delete emPhysicsList;
  delete    decPhysicsList;
 delete radPhysicsList;
  
}


void FlashPhysicsList::AddPackage(const G4String& name)
{
  G4PhysListFactory factory;
  G4VModularPhysicsList* phys =factory.GetReferencePhysList(name);
  G4int i=0;
  const G4VPhysicsConstructor* elem= phys->GetPhysics(i);
  G4VPhysicsConstructor* tmp = const_cast<G4VPhysicsConstructor*> (elem);
  while (elem !=0)
	{
	  RegisterPhysics(tmp);
	  elem= phys->GetPhysics(++i) ;
	  tmp = const_cast<G4VPhysicsConstructor*> (elem);
	}
}




void FlashPhysicsList::ConstructParticle()
{
  decPhysicsList->ConstructParticle();
}

void FlashPhysicsList::ConstructProcess()
{
  // transportation
  //
   AddTransportation();

  // electromagnetic physics list
  //
  emPhysicsList->ConstructProcess();
  em_config.AddModels();

   decPhysicsList->ConstructProcess();
  radPhysicsList->ConstructProcess();
  
  }

void FlashPhysicsList::AddPhysicsList(const G4String& name)
{

  if (verboseLevel>1) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }
  if (name == emName) return;

  /////////////////////////////////////////////////////////////////////////////
  //   ELECTROMAGNETIC MODELS
  /////////////////////////////////////////////////////////////////////////////

  if (name == "standard_opt3") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option3();
    G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option3" << G4endl;
  }
  else if (name == "standard_opt4") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option4();
    G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option4" << G4endl;
    

  } else if (name == "Livermore") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmLivermorePhysics();
    G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmLivermorePhysics" << G4endl;

  } else if (name == "Penelope") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmPenelopePhysics();
    G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmPenelopePhysics" << G4endl;

  }
     
}
void FlashPhysicsList::SetCuts() {//set cuts based on region name defined in detector construction

  SetCutsWithDefault();
  G4Region *region;
  G4String regName;
  G4ProductionCuts *cuts;



  regName = "Phantom_reg";
  region = G4RegionStore::GetInstance()->GetRegion(regName);
  cuts = new G4ProductionCuts;
  cuts->SetProductionCut(0.1 * mm, G4ProductionCuts::GetIndex("gamma"));
  cuts->SetProductionCut(0.1 * mm, G4ProductionCuts::GetIndex("e-"));
  cuts->SetProductionCut(0.1 * mm, G4ProductionCuts::GetIndex("e+"));
  region->SetProductionCuts(cuts);
  
 
}

 

