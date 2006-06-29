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
// $Id: StackingMessenger.cc,v 1.3 2006-06-29 17:24:32 gunter Exp $
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
