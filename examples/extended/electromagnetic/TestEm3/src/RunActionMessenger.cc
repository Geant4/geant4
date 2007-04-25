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
// $Id: RunActionMessenger.cc,v 1.13 2007-04-25 10:07:27 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
:Run(run)
{
  runDir = new G4UIdirectory("/testem/run/");
  runDir->SetGuidance("run commands");
    
  accCmd = new G4UIcommand("/testem/run/acceptance",this);
  accCmd->SetGuidance("Check Edep and RMS of energy deposition for given absorber");
  //
  G4UIparameter* AbsNbPrm = new G4UIparameter("AbsorNb",'i',false);
  AbsNbPrm->SetGuidance("absorber number : from 1 to NbOfAbsor");
  AbsNbPrm->SetParameterRange("AbsorNb>0");
  accCmd->SetParameter(AbsNbPrm);
  //    
  G4UIparameter* edep = new G4UIparameter("Edep",'d',false);
  edep->SetGuidance("mean energy deposition (MeV)");
  edep->SetParameterRange("Edep>=0.");
  accCmd->SetParameter(edep);
  //    
  G4UIparameter* rms = new G4UIparameter("RMS",'d',false);
  rms->SetGuidance("RMS of energy deposition (MeV)");
  rms->SetParameterRange("RMS>=0.");
  accCmd->SetParameter(rms);
  //    
  G4UIparameter* lim = new G4UIparameter("nRMS",'d',false);
  lim->SetGuidance("Limit in number of RMS of energy deposition");
  lim->SetParameterRange("Limit>=0.");
  accCmd->SetParameter(lim);
  //
  limCmd = new G4UIcmdWithABool("/testem/run/limitEdep",this);
  limCmd->SetGuidance("remove energy outside acceptance limit");
  limCmd->AvailableForStates(G4State_PreInit,G4State_Idle);      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMessenger::~RunActionMessenger()
{
  delete accCmd;
  delete runDir;    
  delete limCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{      
  if( command == accCmd )
   { 
     G4int num; 
     G4double edep, rms, lim;
     std::istringstream is(newValue);
     is >> num >> edep >> rms >> lim;
     Run->SetEdepAndRMS(num,edep,rms,lim);
   }

  if( command == limCmd )
   { Run->SetApplyLimit(limCmd->GetNewBoolValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
