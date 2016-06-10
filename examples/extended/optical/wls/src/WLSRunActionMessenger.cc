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
// $Id: WLSRunActionMessenger.cc 70603 2013-06-03 11:23:16Z gcosmo $
//
/// \file optical/wls/src/WLSRunActionMessenger.cc
/// \brief Implementation of the WLSRunActionMessenger class
//
//
#include "globals.hh"
#include "Randomize.hh"

#include "G4UImanager.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

#include "WLSRunAction.hh"
#include "WLSRunActionMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSRunActionMessenger::WLSRunActionMessenger(WLSRunAction* runaction)
  : fRunAction (runaction)
{
  fRndmDir = new G4UIdirectory("/rndm/");
  fRndmDir->SetGuidance("Rndm status control.");
 
  fRndmSaveCmd = new G4UIcmdWithAnInteger("/rndm/save",this);
  fRndmSaveCmd->
           SetGuidance("set frequency to save rndm status on external files.");
  fRndmSaveCmd->SetGuidance("freq = 0 not saved");
  fRndmSaveCmd->SetGuidance("freq > 0 saved on: beginOfRun.rndm");
  fRndmSaveCmd->SetGuidance("freq = 1 saved on:   endOfRun.rndm");
  fRndmSaveCmd->SetGuidance("freq = 2 saved on: endOfEvent.rndm");    
  fRndmSaveCmd->SetParameterName("frequency",false);
  fRndmSaveCmd->SetRange("frequency>=0 && frequency<=2");
  fRndmSaveCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  fRndmReadCmd = new G4UIcmdWithAString("/rndm/read",this);
  fRndmReadCmd->SetGuidance("get rndm status from an external file.");
  fRndmReadCmd->SetParameterName("fileName",true);
  fRndmReadCmd->SetDefaultValue ("beginOfRun.rndm");
  fRndmReadCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSetAutoSeedCmd = new G4UIcmdWithABool("/rndm/autoSeed",this);
  fSetAutoSeedCmd->SetGuidance("Switch on/off time-based random seeds");
  fSetAutoSeedCmd->SetGuidance(" true: run seeds determined by system time");
  fSetAutoSeedCmd->SetGuidance("false: use command 'random/resetEngineFrom'");
  fSetAutoSeedCmd->SetGuidance("Default = false");
  fSetAutoSeedCmd->SetParameterName("autoSeed", false);
  fSetAutoSeedCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSRunActionMessenger::~WLSRunActionMessenger()
{
  delete fRndmDir; delete fRndmSaveCmd;
  delete fRndmReadCmd; delete fSetAutoSeedCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSRunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if (command == fRndmSaveCmd)
      fRunAction->SetRndmFreq(fRndmSaveCmd->GetNewIntValue(newValue));

  if (command == fRndmReadCmd)
  {  G4cout << "\n---> rndm status restored from file: " << newValue << G4endl;
     G4Random::restoreEngineStatus(newValue);
     G4Random::showEngineStatus();
  }

  if(command == fSetAutoSeedCmd)
      fRunAction->SetAutoSeed(fSetAutoSeedCmd->GetNewBoolValue(newValue));
}
