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
#include "G4UIcmdWithAString.hh"
#include "Tst50PrimaryGeneratorMessenger.hh"
#include "Tst50PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst50PrimaryGeneratorMessenger::Tst50PrimaryGeneratorMessenger(Tst50PrimaryGeneratorAction* primary)
:tst50Gun(primary)
{ 
  randomDirectionCmd = new G4UIcmdWithAString("/gun/randomDirection",this);
  randomDirectionCmd -> SetGuidance("Shoot  incident particle with random direction.");
  randomDirectionCmd -> SetGuidance("  Choice : on, off (default)");
  randomDirectionCmd -> SetParameterName("choice",true);
  randomDirectionCmd -> SetDefaultValue("off");
  randomDirectionCmd -> SetCandidates("on off");
  randomDirectionCmd -> AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst50PrimaryGeneratorMessenger::~Tst50PrimaryGeneratorMessenger()
{
  delete  randomDirectionCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst50PrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if( command == randomDirectionCmd )
   { tst50Gun->SetRandomDirection(newValue);} 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

