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

#include "Em3PhysicsList.hh"
#include "Em3PhysicsListMessenger.hh"
 
#include "G4UnitsTable.hh"
#include "Em3PhysListParticles.hh"
#include "Em3PhysListGeneral.hh"
#include "Em3PhysListEmStandard.hh"
#include "Em3PhysListEmModel.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em3PhysicsList::Em3PhysicsList() : G4VModularPhysicsList(),
  emPhysicsListIsRegistered(false)
{   
  currentDefaultCut   = 1.0*mm;
  cutForGamma         = currentDefaultCut;
  cutForElectron      = currentDefaultCut;
  cutForPositron      = currentDefaultCut;

  pMessenger = new Em3PhysicsListMessenger(this);

  SetVerboseLevel(1);

  // Particles
  RegisterPhysics( new Em3PhysListParticles("particles") );

  // General Physics
  RegisterPhysics( new Em3PhysListGeneral("general") );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em3PhysicsList::~Em3PhysicsList()
{
  delete pMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em3PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>1) {
    G4cout << "Em3PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if(name == "standard") {

    if (emPhysicsListIsRegistered) {

      G4cout << "Em3PhysicsList::AddPhysicsList: <" << name << ">" 
             << " cannot be register additionally to existing one"
             << G4endl;
    } else {

      RegisterPhysics( new Em3PhysListEmStandard(name) );
      emPhysicsListIsRegistered = true;
    }

  } else if(name == "model") {

    if (emPhysicsListIsRegistered) {

      G4cout << "Em3PhysicsList::AddPhysicsList: <" << name << ">" 
             << " cannot be register additionally to existing one"
             << G4endl;
    } else {

      RegisterPhysics( new Em3PhysListEmModel(name) );
      emPhysicsListIsRegistered = true;
    }

  } else {

    G4cout << "Em3PhysicsList::AddPhysicsList: <" << name << ">" 
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em3PhysicsList::SetCuts()
{
     
  if (verboseLevel >0){
    G4cout << "Em3PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }  

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");   
    
  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em3PhysicsList::SetCutForGamma(G4double cut)
{
  ResetCuts();
  cutForGamma = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em3PhysicsList::SetCutForElectron(G4double cut)
{
  ResetCuts();
  cutForElectron = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em3PhysicsList::SetCutForPositron(G4double cut)
{
  ResetCuts();
  cutForPositron = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

