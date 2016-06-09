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
// $Id: RunActionMessenger.cc,v 1.8 2004/10/20 14:32:37 maire Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunActionMessenger.hh"

#include "RunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWith3Vector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMessenger::RunActionMessenger(RunAction* run)
:Run(run)
{
  runDir = new G4UIdirectory("/testem/run/");
  runDir->SetGuidance("run control");
    
  accCmd1 = new G4UIcmdWith3Vector("/testem/run/acceptanceL1",this);
  accCmd1->SetGuidance("set Edep and RMS");
  accCmd1->SetGuidance("acceptance values for first layer");
  accCmd1->SetParameterName("edep","rms","limit",true);
  accCmd1->SetRange("edep>0 && edep<1 && rms>0");
  accCmd1->AvailableForStates(G4State_PreInit,G4State_Idle);

  accCmd2 = new G4UIcmdWith3Vector("/testem/run/acceptanceL2",this);
  accCmd2->SetGuidance("set Edep and RMS");
  accCmd2->SetGuidance("acceptance values for 2nd layer");
  accCmd2->SetParameterName("edep","rms","limit",true);
  accCmd2->SetRange("edep>0 && edep<1 && rms>0");
  accCmd2->AvailableForStates(G4State_PreInit,G4State_Idle);      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMessenger::~RunActionMessenger()
{
  delete accCmd1;
  delete accCmd2;
  delete runDir;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{      
  if( command == accCmd1 )
   { Run->SetEdepAndRMS(1,accCmd1->GetNew3VectorValue(newValue));}

  if( command == accCmd2 )
   { Run->SetEdepAndRMS(2,accCmd2->GetNew3VectorValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
