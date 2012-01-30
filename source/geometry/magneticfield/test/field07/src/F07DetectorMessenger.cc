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
//
// $Id: F07DetectorMessenger.cc,v 1.12 2008-09-22 16:41:20 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F07DetectorMessenger.hh"

#include "F07DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F07DetectorMessenger::F07DetectorMessenger(F07DetectorConstruction* myDet)
:myDetector(myDet)
{ 
  F07Dir = new G4UIdirectory("/F07/");
  F07Dir->SetGuidance("UI commands specific to this example.");
  
  detDir = new G4UIdirectory("/F07/det/");
  detDir->SetGuidance("detector control.");
  
  TargMatCmd = new G4UIcmdWithAString("/F07/det/setTargetMate",this);
  TargMatCmd->SetGuidance("Select Material of the Target.");
  TargMatCmd->SetParameterName("choice",false);
  TargMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  ChamMatCmd = new G4UIcmdWithAString("/F07/det/setChamberMate",this);
  ChamMatCmd->SetGuidance("Select Material of the Target.");
  ChamMatCmd->SetParameterName("choice",false);
  ChamMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
  
  FieldCmd = new G4UIcmdWithADoubleAndUnit("/F07/det/setField",this);  
  FieldCmd->SetGuidance("Define magnetic field.");
  FieldCmd->SetGuidance("Magnetic field will be in X direction.");
  FieldCmd->SetParameterName("Bx",false);
  FieldCmd->SetUnitCategory("Magnetic flux density");
  FieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
        
  StepMaxCmd = new G4UIcmdWithADoubleAndUnit("/F07/det/stepMax",this);  
  StepMaxCmd->SetGuidance("Define a step max");
  StepMaxCmd->SetParameterName("stepMax",false);
  StepMaxCmd->SetUnitCategory("Length");
  StepMaxCmd->AvailableForStates(G4State_Idle);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F07DetectorMessenger::~F07DetectorMessenger()
{
  delete TargMatCmd;
  delete ChamMatCmd;
  delete FieldCmd;
  delete StepMaxCmd;  
  delete detDir;
  delete F07Dir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F07DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == TargMatCmd )
   { myDetector->setTargetMaterial(newValue);}
   
  if( command == ChamMatCmd )
   { myDetector->setChamberMaterial(newValue);}  
  
  if( command == FieldCmd )
   { myDetector->SetMagField(FieldCmd->GetNewDoubleValue(newValue));}
      
  if( command == StepMaxCmd )
   { myDetector->SetMaxStep(StepMaxCmd->GetNewDoubleValue(newValue));}   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
