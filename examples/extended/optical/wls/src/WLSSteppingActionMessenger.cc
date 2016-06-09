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
/// \file optical/wls/src/WLSSteppingActionMessenger.cc
/// \brief Implementation of the WLSSteppingActionMessenger class
//
//
//

#include "G4UIdirectory.hh"
#include "WLSSteppingAction.hh"

#include "WLSSteppingActionMessenger.hh"

#include "G4UIcmdWithAnInteger.hh"

WLSSteppingActionMessenger::WLSSteppingActionMessenger(WLSSteppingAction* SA)
  : steppingAction (SA)
{
  steppingDir = new G4UIdirectory("/stepping/");
  steppingDir->SetGuidance("stepping control");

  SetBounceLimitCmd =
                   new G4UIcmdWithAnInteger("/stepping/setBounceLimit", this);
  SetBounceLimitCmd->SetGuidance("Select the maximum number of allowed bounce");
  SetBounceLimitCmd->
              SetGuidance("Set this number to zero if you don't want to limit");
  SetBounceLimitCmd->SetParameterName("limit",false);
  SetBounceLimitCmd->SetRange("limit>=0");
  SetBounceLimitCmd->AvailableForStates(G4State_Idle);
}

WLSSteppingActionMessenger::~WLSSteppingActionMessenger()
{
  delete steppingDir;
  delete SetBounceLimitCmd;
}

void WLSSteppingActionMessenger::SetNewValue(G4UIcommand* command,
                                             G4String newValue)
{
  if ( command == SetBounceLimitCmd ) {

     steppingAction->
               SetBounceLimit(G4UIcmdWithAnInteger::GetNewIntValue(newValue));
  }
}
