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
