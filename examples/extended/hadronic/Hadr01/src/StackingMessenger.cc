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
// $Id: StackingMessenger.cc,v 1.2 2006-06-06 19:48:38 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/////////////////////////////////////////////////////////////////////////
//
// StackingMessenger
//
// Created: 31.05.2006 V.Ivanchenko
//
// Modified:
// 04.06.2006 Adoptation of hadr01 (V.Ivanchenko)
//
////////////////////////////////////////////////////////////////////////
//

#include "StackingMessenger.hh"

#include "StackingAction.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingMessenger::StackingMessenger(StackingAction* stack)
:stackAction(stack)
{
  killCmd = new G4UIcmdWithABool("/testhadr/KillAllSecondaries",this);
  killCmd->SetGuidance("  Choice : true false");
  killCmd->SetParameterName("choice",true);
  killCmd->SetDefaultValue(false);

  kCmd = new G4UIcmdWithAString("/testhadr/Kill", this);
  kCmd->SetGuidance("Kill secondary particles of defined type");
  kCmd->SetParameterName("ch", true);
  kCmd->SetDefaultValue("none");
  kCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingMessenger::~StackingMessenger()
{
  delete killCmd;
  delete kCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StackingMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{     
  if(command == killCmd)
    stackAction->SetKillStatus(killCmd->GetNewBoolValue(newValue));               

  if(command == kCmd)
    stackAction->SetKill(newValue);               
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
