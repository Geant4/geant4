//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "PhysListParticles.hh"
#include "PhysListGeneral.hh"
#include "PhysListEmStandard.hh"
#include "PhysListEmModel.hh"
#include "PhysListEmModelP.hh"
#include "PhysListHadronElastic.hh"
#include "PhysListBinaryCascade.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  G4LossTableManager::Instance();
  defaultCutValue = 1.*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  vertexDetectorCuts = 0;
  muonDetectorCuts   = 0;

  stepMaxProcess  = 0;

  pMessenger = new PhysicsListMessenger(this);

  SetVerboseLevel(1);

   // Particles
  particleList = new PhysListParticles("particles");

  // General Physics
  generalPhysicsList = new PhysListGeneral("general");

  // EM physics
  emName = G4String("standard");
  emPhysicsList = new PhysListEmStandard(emName);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete pMessenger;
  delete generalPhysicsList;
  delete emPhysicsList;
  for(size_t i=0; i<hadronPhys.size(); i++) {
    delete hadronPhys[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  particleList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  generalPhysicsList->ConstructProcess();
  emPhysicsList->ConstructProcess();
  for(size_t i=0; i<hadronPhys.size(); i++) {
    hadronPhys[i]->ConstructProcess();
  }
  AddStepMax();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>-1) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if (name == emName) return;

  if (name == "standard") {

    emName = name;
    delete emPhysicsList;
    emPhysicsList = new PhysListEmStandard(name);

  } else if (name == "model") {

    emName = name;
    delete emPhysicsList;
    emPhysicsList = new PhysListEmModel(name);

  } else if (name == "modelP") {

    emName = name;
    delete emPhysicsList;
    emPhysicsList = new PhysListEmModelP(name);


  } else if (name == "elastic") {

    hadronPhys.push_back( new PhysListHadronElastic(name));

  } else if (name == "binary") {

    hadronPhys.push_back( new PhysListBinaryCascade(name));

  } else {

    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StepMax.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddStepMax()
{
  // Step limitation seen as a process
  stepMaxProcess = new StepMax();

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

  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
G4cout << "world cuts are set" << G4endl;
  if( vertexDetectorCuts ) {
    G4Region* region = (G4RegionStore::GetInstance())->GetRegion("VertexDetector");
    region->SetProductionCuts(vertexDetectorCuts);
  }
G4cout << "Vertex cuts are set" << G4endl;

  if( muonDetectorCuts ) {
    G4Region* region = (G4RegionStore::GetInstance())->GetRegion("MuonDetector");
    region->SetProductionCuts(muonDetectorCuts);
  }
G4cout << "Muon cuts are set" << G4endl;

  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForGamma(G4double cut)
{
  cutForGamma = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForElectron(G4double cut)
{
  cutForElectron = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForPositron(G4double cut)
{
  cutForPositron = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetVertexCut(G4double cut)
{
  if( !vertexDetectorCuts ) vertexDetectorCuts = new G4ProductionCuts();

  vertexDetectorCuts->SetProductionCut(cut, idxG4GammaCut);
  vertexDetectorCuts->SetProductionCut(cut, idxG4ElectronCut);
  vertexDetectorCuts->SetProductionCut(cut, idxG4PositronCut);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetMuonCut(G4double cut)
{
  if( !muonDetectorCuts ) muonDetectorCuts = new G4ProductionCuts();

  muonDetectorCuts->SetProductionCut(cut, idxG4GammaCut);
  muonDetectorCuts->SetProductionCut(cut, idxG4ElectronCut);
  muonDetectorCuts->SetProductionCut(cut, idxG4PositronCut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
