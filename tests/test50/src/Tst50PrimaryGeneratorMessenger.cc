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

