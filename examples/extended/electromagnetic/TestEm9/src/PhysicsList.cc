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
// $Id: PhysicsList.cc,v 1.26 2010-08-12 11:38:20 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   PhysicsList
//
// Author:      V.Ivanchenko 14.10.2002
//
// Modified:
// 17.11.06 Use components from physics_lists subdirectory (V.Ivanchenko)
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "PhysListEmStandard.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4QStoppingPhysics.hh"

#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"
#include "StepMax.hh"

#include "G4EmProcessOptions.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  G4LossTableManager::Instance();
  defaultCutValue = 1.*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;
  cutForVertexDetector = defaultCutValue;
  cutForMuonDetector   = defaultCutValue;

  vertexDetectorCuts = 0;
  muonDetectorCuts   = 0;

  pMessenger = new PhysicsListMessenger(this);
  stepMaxProcess = new StepMax();

  // Initilise flags

  SetVerboseLevel(1);

  helIsRegisted  = false;
  bicIsRegisted  = false;
  gnucIsRegisted = false;
  stopIsRegisted = false;

  // Decay Physics is always defined
  generalPhysicsList = new G4DecayPhysics();

  // EM physics
  emName = G4String("emstandard");
  emPhysicsList = new G4EmStandardPhysics();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete pMessenger;
  delete generalPhysicsList;
  delete emPhysicsList;
  delete stepMaxProcess;
  for(size_t i=0; i<hadronPhys.size(); i++) {
    delete hadronPhys[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  generalPhysicsList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  emPhysicsList->ConstructProcess();
  generalPhysicsList->ConstructProcess();
  for(size_t i=0; i<hadronPhys.size(); i++) {
    hadronPhys[i]->ConstructProcess();
  }
  AddStepMax();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel > 1) 
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;

  if (name == emName) return;

  if (name == "emstandard") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics();
    if (verboseLevel > 0) 
      G4cout << "PhysicsList::Set " << name << " EM physics" << G4endl;

  } else if (name == "emstandard_opt1") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option1();
    if (verboseLevel > 0) 
      G4cout << "PhysicsList::Set " << name << " EM physics" << G4endl;

  } else if (name == "emstandard_opt2") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option2();
    if (verboseLevel > 0) 
      G4cout << "PhysicsList::Set " << name << " EM physics" << G4endl;

  } else if (name == "emstandard_opt3") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option3();
    if (verboseLevel > 0) 
      G4cout << "PhysicsList::Set " << name << " EM physics" << G4endl;

  } else if (name == "emstandard_local") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new PhysListEmStandard();
    if (verboseLevel > 0) 
      G4cout << "PhysicsList::Set " << name << " EM physics" << G4endl;

  } else if (name == "emlivermore") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmLivermorePhysics();

  } else if (name == "empenelope") {
    emName = name;
    delete emPhysicsList;
    emPhysicsList = new G4EmPenelopePhysics();

  } else if (name == "elastic" && !helIsRegisted) {
    hadronPhys.push_back( new G4HadronElasticPhysics());
    helIsRegisted = true;
    if (verboseLevel > 0) 
      G4cout << "PhysicsList::Add hadron elastic physics" << G4endl;

  } else if (name == "binary" && !bicIsRegisted) {
    hadronPhys.push_back(new G4HadronInelasticQBBC());
    hadronPhys.push_back(new G4IonBinaryCascadePhysics());
    bicIsRegisted = true;
    if (verboseLevel > 0) 
      G4cout << "PhysicsList::Add hadron inelastic physics from <QBBC>" << G4endl;

  } else if (name == "gamma_nuc" && !gnucIsRegisted) {
    hadronPhys.push_back(new G4EmExtraPhysics());
    gnucIsRegisted = true;
    if (verboseLevel > 0) 
      G4cout << "PhysicsList::Add gamma- and electro-nuclear physics" << G4endl;

  } else if (name == "stopping" && !stopIsRegisted) {
    hadronPhys.push_back(new G4QStoppingPhysics());
    gnucIsRegisted = true;
    if (verboseLevel > 0) 
      G4cout << "PhysicsList::Add stopping physics" << G4endl;

  } else {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddStepMax()
{
  // Step limitation seen as a process

  theParticleIterator->reset();
  while ((*theParticleIterator)()){
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();

      if (stepMaxProcess->IsApplicable(*particle) && !particle->IsShortLived())
        {
          pmanager ->AddDiscreteProcess(stepMaxProcess);
        }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
  SetCutValue(cutForGamma, "gamma", "DefaultRegionForTheWorld");
  SetCutValue(cutForElectron, "e-", "DefaultRegionForTheWorld");
  SetCutValue(cutForPositron, "e+", "DefaultRegionForTheWorld");
  //  G4cout << "PhysicsList: world cuts are set cutG= " << cutForGamma/mm 
  //	 << " mm    cutE= " << cutForElectron/mm << " mm " << G4endl;

  //G4cout << " cutV= " << cutForVertexDetector 
  //     << " cutM= " << cutForMuonDetector<<G4endl;

  G4Region* region = (G4RegionStore::GetInstance())->GetRegion("VertexDetector");
  vertexDetectorCuts = region->GetProductionCuts();
  SetVertexCut(cutForVertexDetector);
  //  G4cout << "Vertex cuts are set" << G4endl;
 
  region = (G4RegionStore::GetInstance())->GetRegion("MuonDetector");
  muonDetectorCuts = region->GetProductionCuts();
  SetMuonCut(cutForMuonDetector);
  //G4cout << "Muon cuts are set " <<muonRegion << " " << muonDetectorCuts << G4endl;
  
  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForGamma(G4double cut)
{
  cutForGamma = cut;
  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForElectron(G4double cut)
{
  cutForElectron = cut;
  SetParticleCuts(cutForElectron, G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForPositron(G4double cut)
{
  cutForPositron = cut;
  SetParticleCuts(cutForPositron, G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetVertexCut(G4double cut)
{
  cutForVertexDetector = cut;
  
  if( vertexDetectorCuts ) {
    vertexDetectorCuts->SetProductionCut(cut, idxG4GammaCut);
    vertexDetectorCuts->SetProductionCut(cut, idxG4ElectronCut);
    vertexDetectorCuts->SetProductionCut(cut, idxG4PositronCut);
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetMuonCut(G4double cut)
{
  cutForMuonDetector   = cut;

  if( muonDetectorCuts ) {
    muonDetectorCuts->SetProductionCut(cut, idxG4GammaCut);
    muonDetectorCuts->SetProductionCut(cut, idxG4ElectronCut);
    muonDetectorCuts->SetProductionCut(cut, idxG4PositronCut);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

