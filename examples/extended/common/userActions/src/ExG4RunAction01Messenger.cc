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
/// \file ExG4RunAction01Messenger.cc
/// \brief Implementation of the ExG4RunAction01Messenger class

#include "ExG4RunAction01Messenger.hh"
#include "ExG4RunAction01.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ExG4RunAction01Messenger::ExG4RunAction01Messenger(ExG4RunAction01* runAction)
  :  G4UImessenger(),
     fRunAction (runAction),
     fTopDirectory(0),
     fDirectory(0),
     fSetVerboseLevelCmd(0),
     fSetSaveRndmCmd(0),
     fSetReadRndmCmd(0),
     fSetRndmFileNameCmd(0),
     fSetAutoSeedCmd(0)
{
  fTopDirectory = new G4UIdirectory("/ExG4/");
  fTopDirectory->SetGuidance("UI commands of common example classes");
  
  fDirectory = new G4UIdirectory("/ExG4/event/");
  fDirectory->SetGuidance("Event control");
       
  fSetVerboseLevelCmd = new G4UIcmdWithAnInteger("/ExG4/run/verboseLevel",this);
  fSetVerboseLevelCmd->SetGuidance("Set run verbose level ." );
  fSetVerboseLevelCmd->SetParameterName("VerboseLevel",false);
  fSetVerboseLevelCmd->AvailableForStates(G4State_PreInit, G4State_Init);

  fSetSaveRndmCmd = new G4UIcmdWithABool("/ExG4/run/saveRndm",this);
  fSetSaveRndmCmd
    ->SetGuidance("Activate saving random number at beginOfRun and endOfRun .");
  fSetSaveRndmCmd->SetParameterName("SaveRndm",false);
  fSetSaveRndmCmd->AvailableForStates(G4State_PreInit, G4State_Idle); 
         
  fSetReadRndmCmd = new G4UIcmdWithABool("/ExG4/run/readRndm",this);
  fSetReadRndmCmd->SetGuidance("Get rndm status from an external file at beginOfRun.");
  fSetReadRndmCmd->SetParameterName("ReadRndm",false);
  fSetReadRndmCmd->AvailableForStates(G4State_PreInit, G4State_Idle);  

  fSetRndmFileNameCmd = new G4UIcmdWithAString("/ExG4/run/rndmFileName",this);
  fSetRndmFileNameCmd->SetGuidance("Set rndm file name (to read from).");
  fSetRndmFileNameCmd->SetParameterName("RndmFileName",false);
  fSetRndmFileNameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);  

  fSetAutoSeedCmd = new G4UIcmdWithABool("/ExG4/run/autoSeed",this);
  fSetAutoSeedCmd->SetGuidance("Switch on/off time-based random seeds");
  fSetAutoSeedCmd->SetGuidance("  true: run seeds determined by system time");
  fSetAutoSeedCmd->SetGuidance(" false: use command 'random/resetEngineFrom'");
  fSetAutoSeedCmd->SetParameterName("AutoSeed",false);
  fSetAutoSeedCmd->AvailableForStates(G4State_PreInit, G4State_Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ExG4RunAction01Messenger::~ExG4RunAction01Messenger()
{
  delete fTopDirectory;
  delete fDirectory;
  delete fSetVerboseLevelCmd;
  delete fSetSaveRndmCmd; 
  delete fSetReadRndmCmd; 
  delete fSetRndmFileNameCmd; 
  delete fSetAutoSeedCmd; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExG4RunAction01Messenger::SetNewValue(G4UIcommand* command,G4String newValues)
{ 
  if ( command == fSetVerboseLevelCmd ) {
    fRunAction->SetVerboseLevel(fSetVerboseLevelCmd->GetNewIntValue(newValues));
  }
  else if ( command == fSetSaveRndmCmd ) {
    fRunAction->SetSaveRndm(fSetSaveRndmCmd->GetNewBoolValue(newValues));
  }    
  else if ( command == fSetReadRndmCmd ) {
    fRunAction->SetReadRndm(fSetReadRndmCmd->GetNewBoolValue(newValues));
  }   
  else if (command == fSetRndmFileNameCmd) { 
    fRunAction->SetRndmFileName(newValues);
  }   
  else if ( command == fSetAutoSeedCmd ) {
    fRunAction->SetAutoSeed(fSetAutoSeedCmd->GetNewBoolValue(newValues));
  }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
