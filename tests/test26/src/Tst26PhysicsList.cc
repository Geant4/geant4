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
/////////////////////////////////////////////////////////////////////////
//
// test26: Cut per region physics
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
// 14-02-03 Make G4ProductionCuts to be members of the class (V.Ivanchenko)
// 19-02-03 Rename G4ProductionCuts (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst26PhysicsList.hh"
#include "Tst26PhysicsListMessenger.hh"
 
#include "G4UnitsTable.hh"
#include "Tst26PhysListParticles.hh"
#include "Tst26PhysListGeneral.hh"
#include "Tst26PhysListEmStandard.hh"
#include "Tst26PhysListEmModel.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26PhysicsList::Tst26PhysicsList() : G4VModularPhysicsList(),
				       emPhysicsListIsRegistered(false),
				       vertexDetectorCuts(0),
				       muonDetectorCuts(0)
{
  cutForWorld          = 1.0*mm;
  cutForVertexDetector = 0.001*mm;
  cutForMuonDetector   = 10.*mm;

  pMessenger = new Tst26PhysicsListMessenger(this);

  SetVerboseLevel(1);

  // Particles
  RegisterPhysics( new Tst26PhysListParticles("particles") );

  // General Physics
  RegisterPhysics( new Tst26PhysListGeneral("general") );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26PhysicsList::~Tst26PhysicsList()
{
  delete pMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>1) {
    G4cout << "Tst26PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if(name == "standard") {

    if (emPhysicsListIsRegistered) {

      G4cout << "Tst26PhysicsList::AddPhysicsList: <" << name << ">" 
             << " cannot be register additionally to existing one"
             << G4endl;
    } else {

      RegisterPhysics( new Tst26PhysListEmStandard(name) );
      emPhysicsListIsRegistered = true;
    }

  } else if(name == "model") {

    if (emPhysicsListIsRegistered) {

      G4cout << "Tst26PhysicsList::AddPhysicsList: <" << name << ">" 
             << " cannot be register additionally to existing one"
             << G4endl;
    } else {

      RegisterPhysics( new Tst26PhysListEmModel(name) );
      emPhysicsListIsRegistered = true;
    }

  } else {

    G4cout << "Tst26PhysicsList::AddPhysicsList: <" << name << ">" 
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26PhysicsList::SetCuts()
{
     
  if (verboseLevel >0){
    G4cout << "Tst26PhysicsList::SetCuts" << G4endl;
  }  

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForWorld, "gamma");
  SetCutValue(cutForWorld, "e-");
  SetCutValue(cutForWorld, "e+");   

  if( !vertexDetectorCuts ) vertexDetectorCuts = new G4ProductionCuts();
  
  vertexDetectorCuts->SetProductionCut(cutForVertexDetector, idxG4GammaCut);
  vertexDetectorCuts->SetProductionCut(cutForVertexDetector, idxG4ElectronCut);
  vertexDetectorCuts->SetProductionCut(cutForVertexDetector, idxG4PositronCut);

  if( !muonDetectorCuts ) muonDetectorCuts = new G4ProductionCuts();

  muonDetectorCuts->SetProductionCut(cutForMuonDetector, idxG4GammaCut);
  muonDetectorCuts->SetProductionCut(cutForMuonDetector, idxG4ElectronCut);
  muonDetectorCuts->SetProductionCut(cutForMuonDetector, idxG4PositronCut);

  G4Region* region = (G4RegionStore::GetInstance())->GetRegion("VertexDetector");
  region->SetProductionCuts(vertexDetectorCuts);
  region = (G4RegionStore::GetInstance())->GetRegion("MuonDetector");
  region->SetProductionCuts(muonDetectorCuts);

  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26PhysicsList::SetCutForWorld(G4double cut)
{
  ResetCuts();
  cutForWorld = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26PhysicsList::SetCutForVertexDetector(G4double cut)
{
  ResetCuts();
  cutForVertexDetector = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26PhysicsList::SetCutForMuonDetector(G4double cut)
{
  ResetCuts();
  cutForMuonDetector = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

