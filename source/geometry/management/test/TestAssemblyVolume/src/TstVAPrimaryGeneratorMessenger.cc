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
//
// $Id: TstVAPrimaryGeneratorMessenger.cc,v 1.5 2006-06-29 18:34:46 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"

#include "TstVAPrimaryGeneratorMessenger.hh"

#include "TstVAPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"

TstVAPrimaryGeneratorMessenger::TstVAPrimaryGeneratorMessenger(TstVAPrimaryGeneratorAction * myPGA)
:myPGAction(myPGA)
{
  G4bool omitable;

  G4UIdirectory* mydetDir = new G4UIdirectory("/myPGAction/");
  mydetDir->SetGuidance("Primary Generator Action control.");

  standardGun = new G4UIcmdWithoutParameter("/myPGAction/standardGun",this);
  standardGun->SetGuidance("Select standard gun.");

  randomGun = new G4UIcmdWithAString("/myPGAction/randomGun",this);
  randomGun->SetGuidance("Select random gun.");
  randomGun->SetParameterName("random-gun-type", omitable = true);
  randomGun->SetDefaultValue("randomDirection");
  randomGun->SetCandidates
    ("randomDirection randomPosition randomPositionAndDirection");
}

void TstVAPrimaryGeneratorMessenger::SetNewValue
(G4UIcommand * command,G4String newValues)
{
  if( command == standardGun)
  {
    myPGAction->SelectPrimaryGeneratorAction
      (TstVAPrimaryGeneratorAction::standardGun);
  }
  if( command == randomGun)
  {
    if (newValues == "randomDirection") {
      myPGAction->SelectPrimaryGeneratorAction
	(TstVAPrimaryGeneratorAction::randomDirectionGun);
    }
    if (newValues == "randomPosition") {
      myPGAction->SelectPrimaryGeneratorAction
	(TstVAPrimaryGeneratorAction::randomPositionGun);
    }
    if (newValues == "randomPositionAndDirection") {
      myPGAction->SelectPrimaryGeneratorAction
	(TstVAPrimaryGeneratorAction::randomPositionAndDirectionGun);
    }
  }
  return;
}
