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
/// \file SAXSSensitiveDetectorMessenger.cc
/// \brief Implementation of the SAXSSensitiveDetectorMessenger class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SAXSSensitiveDetectorMessenger.hh"
#include "SAXSSensitiveDetector.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SAXSSensitiveDetectorMessenger::SAXSSensitiveDetectorMessenger(SAXSSensitiveDetector* SD):
  G4UImessenger(), fSenDet(SD)
{
  fSenDetDir = new G4UIdirectory("/sd/");
  fSenDetDir->SetGuidance("Sensitive Detector control.");
    
  fUserStopAndKillCmd = new G4UIcmdWithABool("/sd/setSK", this);
  fUserStopAndKillCmd->SetGuidance("bStopAndKill: (1) true, (0) false.");          
  fUserStopAndKillCmd->SetParameterName("bStopAndKill", true);
  fUserStopAndKillCmd->SetDefaultValue(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SAXSSensitiveDetectorMessenger::~SAXSSensitiveDetectorMessenger()
{
  delete fSenDetDir;
  delete fUserStopAndKillCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SAXSSensitiveDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{     
  if (command == fUserStopAndKillCmd) 
    fSenDet->SetVarStopAndKill(fUserStopAndKillCmd->GetNewBoolValue(newValue));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

