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

#include "Tst01PrimaryGeneratorMessenger.hh"

#include "Tst01PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "globals.hh"
#include "G4ios.hh"

Tst01PrimaryGeneratorMessenger::
Tst01PrimaryGeneratorMessenger(Tst01PrimaryGeneratorAction * myPGA)
:fPrimaryGeneratorAction(myPGA)
{
  G4bool omitable;

  G4UIdirectory* mydetDir = new G4UIdirectory("/tsGun/");
  mydetDir->SetGuidance("Primary Generator Action control.");

  standardGun = new G4UIcmdWithoutParameter("/tsGun/standardGun",this);
  standardGun->SetGuidance("Select standard gun.");

  randomGun = new G4UIcmdWithAString("/tsGun/randomGun",this);
  randomGun->SetGuidance("Select random gun.");
  randomGun->SetParameterName("random-gun-type", omitable = true);
  randomGun->SetDefaultValue("randomDirection");
  randomGun->SetCandidates
    ("randomDirection randomPosition randomPositionAndDirection");

  viewerCmd = new G4UIcmdWith3VectorAndUnit("/tsGun/viewerGun",this);
  viewerCmd->SetGuidance("Set starting position of the particle.");
  viewerCmd->SetParameterName("X","Y","Z",true,true);
  viewerCmd->SetDefaultUnit("cm");
  //viewerCmd->SetUnitCategory("Length");
  //viewerCmd->SetUnitCandidates("microm mm cm m km");

  planeCmd = new G4UIcmdWith3VectorAndUnit("/tsGun/planeGun",this);
  planeCmd->SetGuidance("Set starting position of the particle.");
  planeCmd->SetParameterName("X","Y","Z",true,true);
  planeCmd->SetDefaultUnit("cm");
  //planeCmd->SetUnitCategory("Length");
  //planeCmd->SetUnitCandidates("microm mm cm m km");


}

void Tst01PrimaryGeneratorMessenger::SetNewValue( G4UIcommand * command ,
                                                  G4String newValues      )
{
  if( command == standardGun)
  {
    fPrimaryGeneratorAction->SelectPrimaryGeneratorAction
      (Tst01PrimaryGeneratorAction::standardGun);
  }
  if( command == randomGun)
  {
    if (newValues == "randomDirection") 
    {
      fPrimaryGeneratorAction->SelectPrimaryGeneratorAction
	(Tst01PrimaryGeneratorAction::randomDirectionGun);
    }
    if (newValues == "randomPosition") 
    {
      fPrimaryGeneratorAction->SelectPrimaryGeneratorAction
	(Tst01PrimaryGeneratorAction::randomPositionGun);
    }
    if (newValues == "randomPositionAndDirection") 
    {
      fPrimaryGeneratorAction->SelectPrimaryGeneratorAction
	(Tst01PrimaryGeneratorAction::randomPositionAndDirectionGun);
    }
  }
  if(command == viewerCmd)
  {
    fPrimaryGeneratorAction->
    SetGunPosition(viewerCmd->GetNew3VectorValue(newValues)) ;
    fPrimaryGeneratorAction->SelectPrimaryGeneratorAction
	(Tst01PrimaryGeneratorAction::viewerGun);    
  }
  if(command == planeCmd)
  {
    fPrimaryGeneratorAction->
    SetGunPosition(planeCmd->GetNew3VectorValue(newValues)) ;
    fPrimaryGeneratorAction->SelectPrimaryGeneratorAction
	(Tst01PrimaryGeneratorAction::planeGun);    
  }
  return;
}
