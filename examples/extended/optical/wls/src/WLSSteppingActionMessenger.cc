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
// $Id: WLSSteppingActionMessenger.cc 69561 2013-05-08 12:25:56Z gcosmo $
//
/// \file optical/wls/src/WLSSteppingActionMessenger.cc
/// \brief Implementation of the WLSSteppingActionMessenger class
//
//
#include "G4UIdirectory.hh"
#include "WLSSteppingAction.hh"

#include "WLSSteppingActionMessenger.hh"

#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSSteppingActionMessenger::
  WLSSteppingActionMessenger(WLSSteppingAction* steppingaction)
  : fSteppingAction (steppingaction)
{
  fSteppingDir = new G4UIdirectory("/stepping/");
  fSteppingDir->SetGuidance("stepping control");

  fSetBounceLimitCmd =
                   new G4UIcmdWithAnInteger("/stepping/setBounceLimit", this);
  fSetBounceLimitCmd->
                   SetGuidance("Select the maximum number of allowed bounce");
  fSetBounceLimitCmd->
              SetGuidance("Set this number to zero if you don't want to limit");
  fSetBounceLimitCmd->SetParameterName("limit",false);
  fSetBounceLimitCmd->SetRange("limit>=0");
  fSetBounceLimitCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSSteppingActionMessenger::~WLSSteppingActionMessenger()
{
  delete fSteppingDir;
  delete fSetBounceLimitCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSSteppingActionMessenger::SetNewValue(G4UIcommand* command,
                                             G4String newValue)
{
  if ( command == fSetBounceLimitCmd ) {

     fSteppingAction->
               SetBounceLimit(G4UIcmdWithAnInteger::GetNewIntValue(newValue));
  }
}
