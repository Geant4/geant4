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
//

#include "globals.hh"
#include "Randomize.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"

#include "F04RunAction.hh"
#include "F04RunActionMessenger.hh"

F04RunActionMessenger::F04RunActionMessenger(F04RunAction* RA)
  : runAction (RA)
{    
  RndmDir = new G4UIdirectory("/rndm/");
  RndmDir->SetGuidance("Rndm status control.");
  
  RndmSaveCmd = new G4UIcmdWithAnInteger("/rndm/save",this);
  RndmSaveCmd->SetGuidance("set frequency to save rndm status on external files.");
  RndmSaveCmd->SetGuidance("freq = 0 not saved");
  RndmSaveCmd->SetGuidance("freq > 0 saved on: beginOfRun.rndm");
  RndmSaveCmd->SetGuidance("freq = 1 saved on:   endOfRun.rndm");
  RndmSaveCmd->SetGuidance("freq = 2 saved on: endOfEvent.rndm");    
  RndmSaveCmd->SetParameterName("frequency",false);
  RndmSaveCmd->SetRange("frequency>=0 && frequency<=2");
  RndmSaveCmd->AvailableForStates(G4State_PreInit,G4State_Idle); 
         
  RndmReadCmd = new G4UIcmdWithAString("/rndm/read",this);
  RndmReadCmd->SetGuidance("get rndm status from an external file.");
  RndmReadCmd->SetParameterName("fileName",true);
  RndmReadCmd->SetDefaultValue ("beginOfRun.rndm");
  RndmReadCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetAutoSeedCmd = new G4UIcmdWithABool("/rndm/autoSeed",this);
  SetAutoSeedCmd->SetGuidance("Switch on/off time-based random seeds");
  SetAutoSeedCmd->SetGuidance(" true: run seeds determined by system time");
  SetAutoSeedCmd->SetGuidance("false: use command 'random/resetEngineFrom'");
  SetAutoSeedCmd->SetGuidance("Default = false");
  SetAutoSeedCmd->SetParameterName("autoSeed", false);
  SetAutoSeedCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

F04RunActionMessenger::~F04RunActionMessenger()
{ 
  delete RndmDir; delete RndmSaveCmd; delete RndmReadCmd; delete SetAutoSeedCmd;
}

void F04RunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if (command == RndmSaveCmd)
      runAction->SetRndmFreq(RndmSaveCmd->GetNewIntValue(newValue));
		 
  if (command == RndmReadCmd)
  {  G4cout << "\n---> rndm status restored from file: " << newValue << G4endl;
     CLHEP::HepRandom::restoreEngineStatus(newValue);
     CLHEP::HepRandom::showEngineStatus();
  }

  if(command == SetAutoSeedCmd)
      runAction->SetAutoSeed(SetAutoSeedCmd->GetNewBoolValue(newValue));
}
