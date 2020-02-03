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
/// \file exoticphysics/monopole/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * det)
  :G4UImessenger(),
   fDetector(det),
   fDetDir(0),
   fMaterCmd(0),
   fSizeXCmd(0),
   fSizeYZCmd(0),
   fStepSizeCmd(0)
    
{ 
  fDetDir = new G4UIdirectory("/testex/det/");
  fDetDir->SetGuidance("detector construction commands");
      
  fMaterCmd = new G4UIcmdWithAString("/testex/det/setMat",this);
  fMaterCmd->SetGuidance("Select material of the box.");
  fMaterCmd->SetParameterName("choice",false);
  fMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fSizeXCmd = new G4UIcmdWithADoubleAndUnit("/testex/det/setSizeX",this);
  fSizeXCmd->SetGuidance("Set sizeX of the absorber");
  fSizeXCmd->SetParameterName("SizeX",false);
  fSizeXCmd->SetRange("SizeX>0.");
  fSizeXCmd->SetUnitCategory("Length");
  fSizeXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fSizeYZCmd = new G4UIcmdWithADoubleAndUnit("/testex/det/setSizeYZ",this);
  fSizeYZCmd->SetGuidance("Set sizeYZ of the absorber");
  fSizeYZCmd->SetParameterName("SizeYZ",false);
  fSizeYZCmd->SetRange("SizeYZ>0.");
  fSizeYZCmd->SetUnitCategory("Length");
  fSizeYZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fStepSizeCmd = new G4UIcmdWithADoubleAndUnit("/testex/det/setStepSize",this);
  fStepSizeCmd->SetGuidance("Set maxStepSize in the absorber");
  fStepSizeCmd->SetParameterName("StepSize",false);
  fStepSizeCmd->SetRange("StepSize>0.");
  fStepSizeCmd->SetUnitCategory("Length");
  fStepSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fMaterCmd;
  delete fSizeXCmd;
  delete fSizeYZCmd;
  delete fStepSizeCmd;
  delete fDetDir;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fMaterCmd )
   { fDetector->SetMaterial(newValue);}
   
  if( command == fSizeXCmd )
   { fDetector->SetSizeX(fSizeXCmd->GetNewDoubleValue(newValue));}
   
  if( command == fSizeYZCmd )
   { fDetector->SetSizeYZ(fSizeYZCmd->GetNewDoubleValue(newValue));}
      
  if( command == fStepSizeCmd )
   { fDetector->SetMaxStepSize(fStepSizeCmd->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
