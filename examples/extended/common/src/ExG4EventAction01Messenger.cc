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
// $Id$
//
/// \file ExG4EventAction01Messenger.cc
/// \brief Implementation of the ExG4EventAction01Messenger class

#include "ExG4EventAction01Messenger.hh"
#include "ExG4EventAction01.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExG4EventAction01Messenger::ExG4EventAction01Messenger(
                                ExG4EventAction01* eventAction)
  : G4UImessenger(),
    fEventAction(eventAction),
    fTopDirectory(0),
    fDirectory(0),
    fSetVerboseLevelCmd(0),
    fSetSaveRndmCmd(0)
{ 
  fTopDirectory = new G4UIdirectory("/ExG4/");
  fTopDirectory->SetGuidance("UI commands of common example classes");
  
  fDirectory = new G4UIdirectory("/ExG4/event/");
  fDirectory->SetGuidance("Event control");
       
  fSetVerboseLevelCmd 
    = new G4UIcmdWithAnInteger("/ExG4/event/verboseLevel",this);
  fSetVerboseLevelCmd->SetGuidance("Set event verbose level ." );
  fSetVerboseLevelCmd->SetParameterName("VerboseLevel",false);
  fSetVerboseLevelCmd->AvailableForStates(G4State_PreInit, G4State_Init);

  fSetSaveRndmCmd = new G4UIcmdWithABool("/ExG4/event/saveRndm",this);
  fSetSaveRndmCmd
    ->SetGuidance("Activate saving random number at endOfEvent.");
  fSetSaveRndmCmd->SetParameterName("SaveRndm",false);
  fSetSaveRndmCmd->AvailableForStates(G4State_Idle, G4State_Init);     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExG4EventAction01Messenger::~ExG4EventAction01Messenger()
{
  delete fTopDirectory;
  delete fDirectory;
  delete fSetVerboseLevelCmd;
  delete fSetSaveRndmCmd;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExG4EventAction01Messenger::SetNewValue(G4UIcommand* command, 
                                             G4String newValue)
{ 
  if(command == fSetVerboseLevelCmd) {
    fEventAction
      ->SetVerboseLevel(fSetVerboseLevelCmd->GetNewIntValue(newValue));
  }
  else if ( command == fSetSaveRndmCmd ) {
    fEventAction->SetSaveRndm(fSetSaveRndmCmd->GetNewBoolValue(newValue));
  }               
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
