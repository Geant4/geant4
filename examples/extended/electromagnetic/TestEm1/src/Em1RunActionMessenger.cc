// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em1RunActionMessenger.cc,v 1.1 1999-10-11 13:07:48 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em1RunActionMessenger.hh"
#include "Em1RunAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em1RunActionMessenger::Em1RunActionMessenger(Em1RunAction* run)
:Em1Run(run)
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
  RndmSaveCmd->AvailableForStates(PreInit,Idle); 
         
  RndmReadCmd = new G4UIcmdWithAString("/rndm/read",this);
  RndmReadCmd->SetGuidance("get rndm status from an external file.");
  RndmReadCmd->SetParameterName("fileName",true);
  RndmReadCmd->SetDefaultValue ("beginOfRun.rndm");
  RndmReadCmd->AvailableForStates(PreInit,Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em1RunActionMessenger::~Em1RunActionMessenger()
{  
  delete RndmSaveCmd; delete RndmReadCmd; delete RndmDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1RunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if (command == RndmSaveCmd)
      Em1Run->SetRndmFreq(RndmSaveCmd->GetNewIntValue(newValue));
		 
  if (command == RndmReadCmd)
    { G4cout << "\n---> rndm status restored from file: " << newValue << endl;
      HepRandom::restoreEngineStatus(newValue);
      HepRandom::showEngineStatus();
    }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
