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
/// \file electromagnetic/TestEm3/src/RunActionMessenger.cc
/// \brief Implementation of the RunActionMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunActionMessenger.hh"

#include "RunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithABool.hh"

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMessenger::RunActionMessenger(RunAction* run)
:G4UImessenger(),fRunAction(run),
 fRunDir(0),
 fAccCmd(0),
 fLimCmd(0)
{
  fRunDir = new G4UIdirectory("/testem/run/");
  fRunDir->SetGuidance("run commands");
    
  fAccCmd = new G4UIcommand("/testem/run/acceptance",this);
  fAccCmd->SetGuidance("Check Edep and RMS of energy deposition for given absorber");
  //
  G4UIparameter* AbsNbPrm = new G4UIparameter("AbsorNb",'i',false);
  AbsNbPrm->SetGuidance("absorber number : from 1 to NbOfAbsor");
  AbsNbPrm->SetParameterRange("AbsorNb>0");
  fAccCmd->SetParameter(AbsNbPrm);
  //    
  G4UIparameter* edep = new G4UIparameter("Edep",'d',false);
  edep->SetGuidance("mean energy deposition (MeV)");
  edep->SetParameterRange("Edep>=0.");
  fAccCmd->SetParameter(edep);
  //    
  G4UIparameter* rms = new G4UIparameter("RMS",'d',false);
  rms->SetGuidance("RMS of energy deposition (MeV)");
  rms->SetParameterRange("RMS>=0.");
  fAccCmd->SetParameter(rms);
  //    
  G4UIparameter* lim = new G4UIparameter("nRMS",'d',false);
  lim->SetGuidance("Limit in number of RMS of energy deposition");
  lim->SetParameterRange("Limit>=0.");
  fAccCmd->SetParameter(lim);
  //
  fLimCmd = new G4UIcmdWithABool("/testem/run/limitEdep",this);
  fLimCmd->SetGuidance("remove energy outside acceptance limit");
  fLimCmd->AvailableForStates(G4State_PreInit,G4State_Idle);      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMessenger::~RunActionMessenger()
{
  delete fAccCmd;
  delete fRunDir;    
  delete fLimCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{      
  if( command == fAccCmd )
   { 
     G4int num; 
     G4double edep, rms, lim;
     std::istringstream is(newValue);
     is >> num >> edep >> rms >> lim;
     fRunAction->SetEdepAndRMS(num,edep,rms,lim);
   }

  if( command == fLimCmd )
   { fRunAction->SetApplyLimit(fLimCmd->GetNewBoolValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
