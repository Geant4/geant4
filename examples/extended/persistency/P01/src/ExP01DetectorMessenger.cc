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
/// \file persistency/P01/src/ExP01DetectorMessenger.cc
/// \brief Implementation of the ExP01DetectorMessenger class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExP01DetectorMessenger.hh"

#include "ExP01DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExP01DetectorMessenger::ExP01DetectorMessenger(ExP01DetectorConstruction* myDet)
: G4UImessenger(),
  fDetector(myDet),
  fN02Dir(0),
  fDetDir(0),
  fTargMatCmd(0),
  fChamMatCmd(0),    
  fFieldCmd(0)
{ 
  fN02Dir = new G4UIdirectory("/P01/");
  fN02Dir->SetGuidance("UI commands specific to this example.");
  
  fDetDir = new G4UIdirectory("/P01/det/");
  fDetDir->SetGuidance("detector control.");
  
  fTargMatCmd = new G4UIcmdWithAString("/P01/det/setTargetMate",this);
  fTargMatCmd->SetGuidance("Select Material of the Target.");
  fTargMatCmd->SetParameterName("choice",false);
  fTargMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fChamMatCmd = new G4UIcmdWithAString("/P01/det/setChamberMate",this);
  fChamMatCmd->SetGuidance("Select Material of the Target.");
  fChamMatCmd->SetParameterName("choice",false);
  fChamMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
  
  fFieldCmd = new G4UIcmdWithADoubleAndUnit("/P01/det/setField",this);  
  fFieldCmd->SetGuidance("Define magnetic field.");
  fFieldCmd->SetGuidance("Magnetic field will be in X direction.");
  fFieldCmd->SetParameterName("Bx",false);
  fFieldCmd->SetUnitCategory("Magnetic flux density");
  fFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExP01DetectorMessenger::~ExP01DetectorMessenger()
{
  delete fTargMatCmd;
  delete fChamMatCmd;
  delete fFieldCmd;
  delete fDetDir;
  delete fN02Dir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExP01DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fTargMatCmd )
   { fDetector->SetTargetMaterial(newValue);}
   
  if( command == fChamMatCmd )
   { fDetector->SetChamberMaterial(newValue);}  
  
  if( command == fFieldCmd )
   { fDetector->SetMagField(fFieldCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
