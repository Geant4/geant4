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
// This is the *BASIC* version of IORT, a Geant4-based application
//
// Main Authors: G.Russo(a,b), C.Casarino*(c), G.C. Candiano(c), G.A.P. Cirrone(d), F.Romano(d)
// Contributor Authors: S.Guatelli(e)
// Past Authors: G.Arnetta(c), S.E.Mazzaglia(d)
//    
//   (a) Fondazione Istituto San Raffaele G.Giglio, Cefalù, Italy
//   (b) IBFM-CNR , Segrate (Milano), Italy
//   (c) LATO (Laboratorio di Tecnologie Oncologiche), Cefalù, Italy
//   (d) Laboratori Nazionali del Sud of the INFN, Catania, Italy
//   (e) University of Wollongong, Australia
//
//   *Corresponding author, email to carlo.casarino@polooncologicocefalu.it
//////////////////////////////////////////////////////////////////////////////////////////////
//

#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh" 
#include "G4Region.hh"     
#include "G4RegionStore.hh"   
#include "IORTPhysicsList.hh"
#include "IORTPhysicsListMessenger.hh"
#include "IORTStepMax.hh"
#include "G4VPhysicsConstructor.hh"

// Physic lists (contained inside the Geant4 source code, in the 'physicslists folder')
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"  
#include "G4EmPenelopePhysics.hh"   
#include "G4EmExtraPhysics.hh"   
#include "G4ParticleDefinition.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ProcessManager.hh"
#include "globals.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4UnitsTable.hh"
#include "G4DecayPhysics.hh"

IORTPhysicsList::IORTPhysicsList() : G4VModularPhysicsList()
{
  defaultCutValue = 0.1 *mm; 
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");

  DumpCutValuesTable();
  
  stepMaxProcess  = 0;

  pMessenger = new IORTPhysicsListMessenger(this);

  SetVerboseLevel(1);

  // EM physics
  emPhysicsList = new G4EmStandardPhysics_option4(1);
  emName = G4String("emstandard_opt4");
  decPhysicsList = new G4DecayPhysics();
}

IORTPhysicsList::~IORTPhysicsList()
{
  delete pMessenger;
  delete emPhysicsList;
  delete decPhysicsList;
}

void IORTPhysicsList::ConstructParticle()
{
  decPhysicsList->ConstructParticle();
}

void IORTPhysicsList::ConstructProcess()
{
  // transportation
  AddTransportation();

  // electromagnetic physics list
  emPhysicsList->ConstructProcess();
  em_config.AddModels();

  // Set cuts for detector
  SetDetectorCut(defaultCutValue); 
  // step limitation (as a full process)
  //
  AddStepMax();
}

void IORTPhysicsList::AddPhysicsList(const G4String& name)
{

    if (verboseLevel>1) {
	G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
    }
    if (name == emName) return;

  if (name == "standard_opt3") {
	emName = name;
	delete emPhysicsList;
	emPhysicsList = new G4EmStandardPhysics_option3();
	G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
	G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option3" << G4endl;


    } else if (name == "livermore") {   
	emName = name;
	delete emPhysicsList;
	emPhysicsList = new G4EmLivermorePhysics();
	G4RunManager::GetRunManager()-> PhysicsHasBeenModified();
	G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmLivermorePhysics" << G4endl;

    } else if (name == "penelope") {    
	emName = name;
	delete emPhysicsList;
	emPhysicsList = new G4EmPenelopePhysics();
	G4RunManager::GetRunManager()-> PhysicsHasBeenModified();
	G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmPenelopePhysics" << G4endl;}

 else if (name == "standard_opt4") {
	emName = name;
	delete emPhysicsList;
	emPhysicsList = new G4EmStandardPhysics_option4();
	G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
	G4cout << "THE FOLLOWING ELECTROMAGNETIC PHYSICS LIST HAS BEEN ACTIVATED: G4EmStandardPhysics_option4" << G4endl;
        }

    else {

	G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
	    << " is not defined"
	    << G4endl;
    }
}

void IORTPhysicsList::AddStepMax()
{
  // Step limitation seen as a process
  stepMaxProcess = new IORTStepMax();

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)()){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();

    if (stepMaxProcess->IsApplicable(*particle) && pmanager)
      {
	pmanager ->AddDiscreteProcess(stepMaxProcess);
      }
  }
}

void IORTPhysicsList::SetCutForGamma(G4double cut)
{
  cutForGamma = cut;
  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

void IORTPhysicsList::SetCutForElectron(G4double cut)
{
  cutForElectron = cut;
  SetParticleCuts(cutForElectron, G4Electron::Electron());
}

void IORTPhysicsList::SetCutForPositron(G4double cut)
{
  cutForPositron = cut;
  SetParticleCuts(cutForPositron, G4Positron::Positron());
}

void IORTPhysicsList::SetDetectorCut(G4double cut)
{
  G4String regionName = "DetectorLog";
  G4Region* region = G4RegionStore::GetInstance()->GetRegion(regionName);

  G4ProductionCuts* cuts = new G4ProductionCuts ;
  cuts -> SetProductionCut(cut,"gamma");
  cuts -> SetProductionCut(cut,"e-");
  cuts -> SetProductionCut(cut,"e+");
  region -> SetProductionCuts(cuts);
}

