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

#include "test31PhysicsList.hh"
#include "test31PhysicsListMessenger.hh"
 
#include "G4UnitsTable.hh"
#include "test31Particles.hh"
#include "test31GeneralPhysics.hh"
#include "test31StandardEM.hh"
#include "test31ModelEM.hh"
#include "test31LowEPhysicsList.hh"
#include "test31HadronElastic.hh"
#include "test31PreCompound.hh"
#include "test31LEparametrised.hh"
#include "test31CHIPS.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

test31PhysicsList::test31PhysicsList() : G4VModularPhysicsList(),
  emPhysicsListIsRegistered(false),
  hadPhysicsListIsRegistered(false)
{   
  currentDefaultCut   = 1.0*mm;
  cutForGamma         = currentDefaultCut;
  cutForElectron      = currentDefaultCut;
  cutForPositron      = currentDefaultCut;
  theMaxStep          = DBL_MAX;

  nuclStop = true;
  barkas   = true;
  table    = "";
  verbose  = 1;

  pMessenger = new test31PhysicsListMessenger(this);

  SetVerboseLevel(verbose);

  // Particles
  RegisterPhysics( new test31Particles("particles") );

  // General Physics
  RegisterPhysics( new test31GeneralPhysics("general", &theStepCut) );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

test31PhysicsList::~test31PhysicsList()
{
  delete pMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void test31PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>1) {
    G4cout << "test31PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if("standard" == name) {

    if (emPhysicsListIsRegistered) {

      G4cout << "test31PhysicsList::AddPhysicsList: <" << name << ">" 
             << " cannot be register additionally to existing one"
             << G4endl;
    } else {

      RegisterPhysics( new test31StandardEM(name) );
      emPhysicsListIsRegistered = true;
    }

  } else if("model" == name) {

    if (emPhysicsListIsRegistered) {

      G4cout << "test31PhysicsList::AddPhysicsList: <" << name << ">" 
             << " cannot be register additionally to existing one"
             << G4endl;
    } else {

      RegisterPhysics( new test31ModelEM(name) );
      emPhysicsListIsRegistered = true;
    }

  } else if("lowenergy" == name) {

    if (emPhysicsListIsRegistered) {

      G4cout << "test31PhysicsList::AddPhysicsList: <" << name << ">" 
             << " cannot be register additionally to existing one"
             << G4endl;
    } else {

      RegisterPhysics( new test31LowEPhysicsList(name, verbose, nuclStop, barkas, table) );
      emPhysicsListIsRegistered = true;
    }

  } else if("LEparametrised" == name) {

    if (hadPhysicsListIsRegistered) {

      G4cout << "test31PhysicsList::AddPhysicsList: <" << name << ">" 
             << " cannot be register additionally to existing one"
             << G4endl;
    } else {

      RegisterPhysics( new test31LEparametrised(name) );
      hadPhysicsListIsRegistered = true;
    }

  } else if("PreCompound" == name) {

    if (hadPhysicsListIsRegistered) {

      G4cout << "test31PhysicsList::AddPhysicsList: <" << name << ">" 
             << " cannot be register additionally to existing one"
             << G4endl;
    } else {

      RegisterPhysics( new test31PreCompound(name) );
      hadPhysicsListIsRegistered = true;
    }

  } else if("CHIPS" == name) {

    if (hadPhysicsListIsRegistered) {

      G4cout << "test31PhysicsList::AddPhysicsList: <" << name << ">" 
             << " cannot be register additionally to existing one"
             << G4endl;
    } else {

      RegisterPhysics( new test31CHIPS(name) );
      hadPhysicsListIsRegistered = true;
    }
  } else if("elastic" == name) {

    RegisterPhysics( new test31HadronElastic(name) );

  } else {

    G4cout << "test31PhysicsList::AddPhysicsList: <" << name << ">" 
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void test31PhysicsList::SetCuts()
{
     
  if (verboseLevel >0){
    G4cout << "test31PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }  

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");   
  
  // set cut values for proton and anti_proton before all other hadrons
  // because some processes for hadrons need cut values for proton/anti_proton
  SetCutValue(currentDefaultCut, "proton");
  SetCutValue(currentDefaultCut, "anti_proton");
     
  SetCutValueForOthers(currentDefaultCut);
  
  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void test31PhysicsList::SetCutForGamma(G4double cut)
{
  ResetCuts();
  cutForGamma = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void test31PhysicsList::SetCutForElectron(G4double cut)
{
  ResetCuts();
  cutForElectron = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void test31PhysicsList::SetCutForPositron(G4double cut)
{
  ResetCuts();
  cutForPositron = cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

