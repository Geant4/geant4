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
// $Id: ActionInitializer.cc 66241 2012-12-13 18:34:42Z gunter $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ActionInitializer.hh"

#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "StackingAction.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitializer::ActionInitializer() : 
  G4VUserActionInitialization()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitializer::Build() const 
{
  //G4cout << "Sono in ActionInitializer::Build()" << G4endl;

  const DetectorConstruction* detector = 
    static_cast<const DetectorConstruction*>
    (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  //  G4cout << "RunManagerType: " << G4RunManager::GetRunManager()->GetRunManagerType() << G4endl;
  
  // primary generator
  PrimaryGeneratorAction* primary = new PrimaryGeneratorAction(detector);
  SetUserAction(primary);
    

  //Thread-local RunAction: same class, but code controlled by IsMaster()
  SetUserAction(new RunAction());
  SetUserAction(new EventAction());
  SetUserAction(new SteppingAction());
  SetUserAction(new StackingAction());

 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitializer::BuildForMaster() const
{ 
  //Thread-local RunAction: same class, but code controlled by IsMaster()
  RunAction* runAction = new RunAction();
  SetUserAction(runAction);
 
 //G4cout << "Sono in ActionInitializer::BuildForMaster()" << G4endl;
}

