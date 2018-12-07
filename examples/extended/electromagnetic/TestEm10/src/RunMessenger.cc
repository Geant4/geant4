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
/// \file electromagnetic/TestEm10/src/RunMessenger.cc
/// \brief Implementation of the RunMessenger class
//
//
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "RunMessenger.hh"

#include "RunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunMessenger::RunMessenger(RunAction* runAction)
 : G4UImessenger(), 
   fRunAction(runAction), 
   fRndmDir(0), fRndmSaveCmd(0), fRndmReadCmd(0)
{
  fRndmDir = new G4UIdirectory("/rndm/");
  fRndmDir->SetGuidance("Rndm status control.");
  
  fRndmSaveCmd = new G4UIcmdWithAnInteger("/rndm/save",this);
  fRndmSaveCmd->SetGuidance("set frequency to save rndm status on external files.");
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunMessenger::~RunMessenger()
{
  delete fRndmSaveCmd; 
  delete fRndmReadCmd; 
  delete fRndmDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  if (command == fRndmSaveCmd) {
    fRunAction->SetRndmFreq(fRndmSaveCmd->GetNewIntValue(newValues));
  }
                 
  if (command == fRndmReadCmd) {
    G4cout << "\n---> rndm status restored from file: " << newValues << G4endl;
    CLHEP::HepRandom::restoreEngineStatus(newValues);
    CLHEP::HepRandom::showEngineStatus();
  }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

   
