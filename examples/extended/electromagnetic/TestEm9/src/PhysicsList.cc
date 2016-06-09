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
/// \file electromagnetic/TestEm9/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   PhysicsList
//
// Author:      V.Ivanchenko 14.10.2002
//
// Modified:
// 17.11.06 Use components from physics_lists subdirectory (V.Ivanchenko)
// 24.10.12 Migrated to new stopping and ion physics (A.Ribon)
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
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4IonPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"

#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"
#include "StepMax.hh"

#include "G4EmProcessOptions.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  G4LossTableManager::Instance();
  defaultCutValue  = 1.*mm;
  fCutForGamma     = defaultCutValue;
  fCutForElectron  = defaultCutValue;
  fCutForPositron  = defaultCutValue;
  fCutForVertexDetector = defaultCutValue;
  fCutForMuonDetector   = defaultCutValue;

  fVertexDetectorCuts = 0;
  fMuonDetectorCuts   = 0;

  fMessenger = new PhysicsListMessenger(this);
  fStepMaxProcess = new StepMax();

  // Initilise flags

  SetVerboseLevel(1);

  fHelIsRegisted  = false;
  fBicIsRegisted  = false;
  fGnucIsRegisted = false;
  fStopIsRegisted = false;

  // EM physics
  fEmName = G4String("emstandard");
  fEmPhysicsList = new G4EmStandardPhysics();

  // Decay Physics is always defined
  fDecayPhysicsList = new G4DecayPhysics();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete fMessenger;
  delete fDecayPhysicsList;
  delete fEmPhysicsList;
  delete fStepMaxProcess;
  for(size_t i=0; i<fHadronPhys.size(); i++) {
    delete fHadronPhys[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  fDecayPhysicsList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  fEmPhysicsList->ConstructProcess();
  fDecayPhysicsList->ConstructProcess();
  for(size_t i=0; i<fHadronPhys.size(); ++i) {
    fHadronPhys[i]->ConstructProcess();
  }
  AddStepMax();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel > 1) 
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;

  if (name == fEmName) return;

  if (name == "emstandard") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics();
    if (verboseLevel > 0) 
      G4cout << "PhysicsList::Set " << name << " EM physics" << G4endl;

  } else if (name == "emstandard_opt1") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option1();
    if (verboseLevel > 0) 
      G4cout << "PhysicsList::Set " << name << " EM physics" << G4endl;

  } else if (name == "emstandard_opt2") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option2();
    if (verboseLevel > 0) 
      G4cout << "PhysicsList::Set " << name << " EM physics" << G4endl;

  } else if (name == "emstandard_opt3") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option3();
    if (verboseLevel > 0) 
      G4cout << "PhysicsList::Set " << name << " EM physics" << G4endl;

  } else if (name == "emstandard_opt4") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option4();
    if (verboseLevel > 0) 
      G4cout << "PhysicsList::Set " << name << " EM physics" << G4endl;

  } else if (name == "emstandard_local") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new PhysListEmStandard();
    if (verboseLevel > 0) 
      G4cout << "PhysicsList::Set " << name << " EM physics" << G4endl;

  } else if (name == "emlivermore") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLivermorePhysics();

  } else if (name == "empenelope") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmPenelopePhysics();

  } else if (name == "elastic" && !fHelIsRegisted) {
    fHadronPhys.push_back( new G4HadronElasticPhysics());
    fHelIsRegisted = true;
    if (verboseLevel > 0) 
      G4cout << "PhysicsList::Add hadron elastic physics" << G4endl;

  } else if (name == "binary" && !fBicIsRegisted) {
    fHadronPhys.push_back(new G4HadronInelasticQBBC());
    fHadronPhys.push_back(new G4IonPhysics());
    fBicIsRegisted = true;
    if (verboseLevel > 0) 
      G4cout << "PhysicsList::Add hadron inelastic physics from <QBBC>" << G4endl;

  } else if (name == "gamma_nuc" && !fGnucIsRegisted) {
    fHadronPhys.push_back(new G4EmExtraPhysics());
    fGnucIsRegisted = true;
    if (verboseLevel > 0) 
      G4cout << "PhysicsList::Add gamma- and electro-nuclear physics" << G4endl;

  } else if (name == "stopping" && !fStopIsRegisted) {
    fHadronPhys.push_back(new G4StoppingPhysics());
    fStopIsRegisted = true;
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

    if (fStepMaxProcess->IsApplicable(*particle) && !particle->IsShortLived())
      {
        pmanager ->AddDiscreteProcess(fStepMaxProcess);
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
  SetCutValue(fCutForGamma, "gamma", "DefaultRegionForTheWorld");
  SetCutValue(fCutForElectron, "e-", "DefaultRegionForTheWorld");
  SetCutValue(fCutForElectron, "e+", "DefaultRegionForTheWorld");
  SetCutValue(fCutForProton, "proton", "DefaultRegionForTheWorld");
  //  G4cout << "PhysicsList: world cuts are set cutG= " << cutForGamma/mm 
  //         << " mm    cutE= " << fCutForElectron/mm << " mm " << G4endl;

  //G4cout << " cutV= " << fCutForVertexDetector 
  //     << " cutM= " << fCutForMuonDetector<<G4endl;

  G4Region* region = (G4RegionStore::GetInstance())->GetRegion("VertexDetector");
  fVertexDetectorCuts = region->GetProductionCuts();
  SetVertexCut(fCutForVertexDetector);
  //  G4cout << "Vertex cuts are set" << G4endl;
 
  region = (G4RegionStore::GetInstance())->GetRegion("MuonDetector");
  fMuonDetectorCuts = region->GetProductionCuts();
  SetMuonCut(fCutForMuonDetector);
  //G4cout << "Muon cuts are set " <<muonRegion << " " << muonDetectorCuts << G4endl;
  
  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForGamma(G4double cut)
{
  fCutForGamma = cut;
  SetParticleCuts(fCutForGamma, G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForElectron(G4double cut)
{
  fCutForElectron = cut;
  SetParticleCuts(fCutForElectron, G4Electron::Electron());
  SetParticleCuts(fCutForElectron, G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForProton(G4double cut)
{
  fCutForProton = cut;
  SetParticleCuts(fCutForProton, G4Proton::Proton());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetVertexCut(G4double cut)
{
  fCutForVertexDetector = cut;
  
  if( fVertexDetectorCuts ) {
    fVertexDetectorCuts->SetProductionCut(cut, idxG4GammaCut);
    fVertexDetectorCuts->SetProductionCut(cut, idxG4ElectronCut);
    fVertexDetectorCuts->SetProductionCut(cut, idxG4PositronCut);
    fVertexDetectorCuts->SetProductionCut(cut, idxG4ProtonCut);
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetMuonCut(G4double cut)
{
  fCutForMuonDetector = cut;

  if( fMuonDetectorCuts ) {
    fMuonDetectorCuts->SetProductionCut(cut, idxG4GammaCut);
    fMuonDetectorCuts->SetProductionCut(cut, idxG4ElectronCut);
    fMuonDetectorCuts->SetProductionCut(cut, idxG4PositronCut);
    fMuonDetectorCuts->SetProductionCut(cut, idxG4ProtonCut);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

