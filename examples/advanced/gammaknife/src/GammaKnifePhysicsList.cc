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

#include "G4RunManager.hh"
#include "globals.hh"
#include "G4ProcessManager.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4PhysListFactory.hh"
#include "GammaKnifePhysicsList.hh"
#include "GammaKnifePhysicsListMessenger.hh"
#include "GammaKnifeParticles.hh"
#include "G4SystemOfUnits.hh"


// Physic lists (contained inside the Geant4 source code, in the 'physicslists folder')
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmExtraPhysics.hh"

#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

#include "G4LossTableManager.hh"

GammaKnifePhysicsList::GammaKnifePhysicsList(): G4VModularPhysicsList()
{
  G4LossTableManager::Instance();
  defaultCutValue = 1.*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  radioactiveDecayIsRegistered = false;


  messenger = new GammaKnifePhysicsListMessenger(this);
  SetVerboseLevel(1);

  // EM physics
  emPhysicsList = new G4EmStandardPhysics_option3(1);
  emName = G4String("emstandard_opt3");

  // Decay physics and all particles
  decPhysicsList = new G4DecayPhysics();
}

GammaKnifePhysicsList::~GammaKnifePhysicsList()
{ 
  delete messenger;
  delete emPhysicsList;
  delete decPhysicsList;
  for(size_t i=0; i<hadronPhys.size(); i++) 
    delete hadronPhys[i];
}

void GammaKnifePhysicsList::ConstructParticle()
{
  decPhysicsList->ConstructParticle();
  emPhysicsList->ConstructParticle();
}

void GammaKnifePhysicsList::ConstructProcess()
{
  // transportation
  AddTransportation();

  // electromagnetic physics list
  emPhysicsList->ConstructProcess();
  
  // decay physics list
  decPhysicsList->ConstructProcess();

  // hadronic physics lists
  for(size_t i=0; i<hadronPhys.size(); i++) {
    hadronPhys[i] -> ConstructProcess();
  }

}


void GammaKnifePhysicsList::AddPhysicsList(const G4String& name)
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
	G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
	G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option3" << G4endl;

    } else if (name == "LowE_Livermore") {
	emName = name;
	delete emPhysicsList;
	emPhysicsList = new G4EmLivermorePhysics();
	G4RunManager::GetRunManager()-> PhysicsHasBeenModified();
	G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmLivermorePhysics" << G4endl;

    } else if (name == "LowE_Penelope") {
	emName = name;
	delete emPhysicsList;
	emPhysicsList = new G4EmPenelopePhysics();
	G4RunManager::GetRunManager()-> PhysicsHasBeenModified();
	G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmLivermorePhysics" << G4endl;

	/////////////////////////////////////////////////////////////////////////////
	//   HADRONIC MODELS
	/////////////////////////////////////////////////////////////////////////////


    }  else if (name == "decay")
    {
	hadronPhys.push_back(new G4DecayPhysics());
	//radioactiveDecayIsRegistered = true;
    }
    else if (name == "radioactive_decay" && !radioactiveDecayIsRegistered )
    {
	hadronPhys.push_back(new G4RadioactiveDecayPhysics());
	radioactiveDecayIsRegistered = true;
    }
    else {

	G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
	    << " is not defined"
	    << G4endl;
    }

 }


/////////////////////////////////////////////////////////////////////////////
void GammaKnifePhysicsList::SetCuts()
{  

  if (verboseLevel >0){
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");

 
  if (verboseLevel>0) DumpCutValuesTable();
}

/////////////////////////////////////////////////////////////////////////////
void GammaKnifePhysicsList::SetCutForGamma(G4double cut)
{
  cutForGamma = cut;
  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

/////////////////////////////////////////////////////////////////////////////
void GammaKnifePhysicsList::SetCutForElectron(G4double cut)
{
  cutForElectron = cut;
  SetParticleCuts(cutForElectron, G4Electron::Electron());
}

/////////////////////////////////////////////////////////////////////////////
void GammaKnifePhysicsList::SetCutForPositron(G4double cut)
{
  cutForPositron = cut;
  SetParticleCuts(cutForPositron, G4Positron::Positron());
}

