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
// $Id: DetectorMessenger.cc,v 1.3 2007-01-18 09:07:20 hbu Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:Detector(Det)
{
  testemDir = new G4UIdirectory("/testem/");
  testemDir->SetGuidance(" detector control.");

  detDir = new G4UIdirectory("/testem/det/");
  detDir->SetGuidance("detector construction commands");

  MaterCmd = new G4UIcmdWithAString("/testem/det/setMat",this);
  MaterCmd->SetGuidance("Select material of the box.");
  MaterCmd->SetParameterName("choice",false);
  MaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SizeCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSize",this);
  SizeCmd->SetGuidance("Set size of the box");
  SizeCmd->SetParameterName("Size",false);
  SizeCmd->SetRange("Size>0.");
  SizeCmd->SetUnitCategory("Length");
  SizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  MagFieldCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setField",this);
  MagFieldCmd->SetGuidance("Define magnetic field.");
  MagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  MagFieldCmd->SetParameterName("Bz",false);
  MagFieldCmd->SetUnitCategory("Magnetic flux density");
  MagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);

  MaxStepCmd = new G4UIcmdWithADoubleAndUnit("/testem/tracking/stepMax",this);
  MaxStepCmd->SetGuidance("Set max allowed step size");
  MaxStepCmd->SetParameterName("Size",false);
  MaxStepCmd->SetRange("Size>0.");
  MaxStepCmd->SetUnitCategory("Length");
  MaxStepCmd->AvailableForStates(G4State_Idle);

  trackdir = new G4UIdirectory("/testem/tracking/");
  trackdir->SetGuidance("step length");
  MaxStepLength = new G4UIcmdWithADoubleAndUnit("/testem/tracking/setMaxStepLength",this);
  MaxStepLength->SetGuidance("Set the maximum length of tracking step");
  MaxStepLength->SetGuidance("when integrating magnetic field line.");
  MaxStepLength->SetParameterName("Size",false);
  MaxStepLength->SetRange("Size>0.");
  MaxStepLength->SetUnitCategory("Length");
  MaxStepLength->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete MaterCmd;
  delete SizeCmd;
  delete MagFieldCmd;
  delete UpdateCmd;
  delete MaxStepCmd;
  delete MaxStepLength;
  delete detDir;
  delete testemDir;
  delete trackdir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == MaterCmd )
   { Detector->SetMaterial(newValue);}

  if( command == SizeCmd )
   { Detector->SetSize(SizeCmd->GetNewDoubleValue(newValue));}

  if( command == MagFieldCmd )
   { Detector->SetMagField(MagFieldCmd->GetNewDoubleValue(newValue));}

  if( command == UpdateCmd )
   { Detector->UpdateGeometry(); }

  if( command == MaxStepCmd )
   { Detector->SetMaxStepSize(MaxStepCmd->GetNewDoubleValue(newValue));}

  if( command == MaxStepLength )
   { Detector->SetMaxStepLength(MaxStepLength->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
