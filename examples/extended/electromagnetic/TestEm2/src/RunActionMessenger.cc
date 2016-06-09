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
// $Id: RunActionMessenger.cc,v 1.1 2004/06/18 15:43:41 maire Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunActionMessenger.hh"

#include "RunAction.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMessenger::RunActionMessenger(RunAction* run)
:Run(run)
{    
  accCmd = new G4UIcmdWith3Vector("/testem/run/acceptance",this);
  accCmd->SetGuidance("set Edep and RMS");
  accCmd->SetGuidance("acceptance values for first layer");
  accCmd->SetParameterName("edep","rms","limit",true);
  accCmd->SetRange("edep>0 && edep<1 && rms>0");
  accCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  factoryCmd = new G4UIcmdWithAString("/testem/histo/fileName",this);
  factoryCmd->SetGuidance("set name for the histograms file");

  fileCmd = new G4UIcmdWithAString("/testem/histo/fileType",this);
  fileCmd->SetGuidance("set type (hbook, root, XML) for the histograms file");         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMessenger::~RunActionMessenger()
{
  delete accCmd;
  delete factoryCmd;
  delete fileCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{      
  if (command == accCmd )
   { Run->SetEdepAndRMS(accCmd->GetNew3VectorValue(newValue));}
   
  if (command == factoryCmd)
   { Run->SetHistoName(newValue);}

  if (command == fileCmd)
   { Run->SetHistoType(newValue);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
