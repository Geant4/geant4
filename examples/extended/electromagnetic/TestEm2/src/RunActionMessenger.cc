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
/// \file electromagnetic/TestEm2/src/RunActionMessenger.cc
/// \brief Implementation of the RunActionMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunActionMessenger.hh"

#include "RunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMessenger::RunActionMessenger(RunAction* run)
:G4UImessenger(),fRun(run),
 fRunDir(0),
 fAccCmd(0),
 fVerbCmd(0), 
 fHistoDir(0),        
 fHFileCmd(0)
{
  fRunDir = new G4UIdirectory("/testem/run/");
  fRunDir->SetGuidance("run control");
      
  fAccCmd = new G4UIcmdWith3Vector("/testem/run/acceptance",this);
  fAccCmd->SetGuidance("set Edep and RMS");
  fAccCmd->SetGuidance("acceptance values for first layer");
  fAccCmd->SetParameterName("edep","rms","limit",true);
  fAccCmd->SetRange("edep>0 && edep<1 && rms>0");
  fAccCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fVerbCmd = new G4UIcmdWithAnInteger("/testem/run/verbose",this);
  fVerbCmd->SetGuidance("set verbose level for runAction");
  fVerbCmd->SetParameterName("verbose",false);
    
  fHistoDir = new G4UIdirectory("/testem/histo/");
  fHistoDir->SetGuidance("histograms control");
  
  fHFileCmd = new G4UIcmdWithAString("/testem/histo/setFileName",this);
  fHFileCmd->SetGuidance("set name for the histograms file");    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMessenger::~RunActionMessenger()
{
  delete fVerbCmd;
  delete fAccCmd;
  delete fRunDir;
  delete fHFileCmd;
  delete fHistoDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{      
  if (command == fAccCmd )
   { fRun->SetEdepAndRMS(fAccCmd->GetNew3VectorValue(newValue));}
   
  if (command == fVerbCmd )
   { fRun->SetVerbose(fVerbCmd->GetNewIntValue(newValue));}
      
  if (command == fHFileCmd)
   { fRun->SetHistoName(newValue);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
