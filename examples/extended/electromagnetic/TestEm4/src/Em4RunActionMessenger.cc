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
//
// $Id: Em4RunActionMessenger.cc,v 1.6 2001-10-17 14:04:16 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em4RunActionMessenger.hh"
#include "Em4RunAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em4RunActionMessenger::Em4RunActionMessenger(Em4RunAction* run)
:Em4Run(run)
{   
  RndmDir = new G4UIdirectory("/rndm/");
  RndmDir->SetGuidance("Rndm status control.");
  
  RndmSaveCmd = new G4UIcmdWithAnInteger("/rndm/save",this);
  RndmSaveCmd->SetGuidance("set frequency to save rndm on external files.");
  RndmSaveCmd->SetGuidance("freq = 0 not saved");
  RndmSaveCmd->SetGuidance("freq > 0 saved on: beginOfRun.rndm");
  RndmSaveCmd->SetGuidance("freq > 0 saved on:   endOfRun.rndm");
  RndmSaveCmd->SetGuidance("freq = 2 saved on: beginOfEvent.rndm");    
  RndmSaveCmd->SetParameterName("frequency",false);
  RndmSaveCmd->SetRange("frequency>=0 && frequency<=2");
  RndmSaveCmd->AvailableForStates(PreInit,Idle); 
         
  RndmReadCmd = new G4UIcmdWithAString("/rndm/read",this);
  RndmReadCmd->SetGuidance("get rndm status from an external file.");
  RndmReadCmd->SetParameterName("fileName",true);
  RndmReadCmd->SetDefaultValue ("beginOfRun.rndm");
  RndmReadCmd->AvailableForStates(PreInit,Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em4RunActionMessenger::~Em4RunActionMessenger()
{  
  delete RndmSaveCmd; delete RndmReadCmd; delete RndmDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em4RunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if (command == RndmSaveCmd)
      Em4Run->SetRndmFreq(RndmSaveCmd->GetNewIntValue(newValue));
		 
  if (command == RndmReadCmd)
    { G4cout << "\n---> rndm status restored from file: " << newValue << G4endl;
      HepRandom::restoreEngineStatus(newValue);
      HepRandom::showEngineStatus();
    }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
