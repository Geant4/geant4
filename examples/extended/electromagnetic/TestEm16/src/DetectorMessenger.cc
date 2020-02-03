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
/// \file electromagnetic/TestEm16/src/DetectorMessenger.cc
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
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(),fDetector(Det),
 fTestemDir(0),
 fDetDir(0),
 fTrackdir(0),
 fMaterCmd(0),
 fSizeCmd(0),
 fMaxStepCmd(0),
 fMaxStepLength(0)
{
  fTestemDir = new G4UIdirectory("/testem/");
  fTestemDir->SetGuidance(" detector control.");

  fDetDir = new G4UIdirectory("/testem/det/");
  fDetDir->SetGuidance("detector construction commands");

  fMaterCmd = new G4UIcmdWithAString("/testem/det/setMat",this);
  fMaterCmd->SetGuidance("Select material of the box.");
  fMaterCmd->SetParameterName("choice",false);
  fMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSizeCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSize",this);
  fSizeCmd->SetGuidance("Set size of the box");
  fSizeCmd->SetParameterName("Size",false);
  fSizeCmd->SetRange("Size>0.");
  fSizeCmd->SetUnitCategory("Length");
  fSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fMaxStepCmd = new G4UIcmdWithADoubleAndUnit("/testem/tracking/stepMax",this);
  fMaxStepCmd->SetGuidance("Set max allowed step size");
  fMaxStepCmd->SetParameterName("Size",false);
  fMaxStepCmd->SetRange("Size>0.");
  fMaxStepCmd->SetUnitCategory("Length");
  fMaxStepCmd->AvailableForStates(G4State_Idle);

  fTrackdir = new G4UIdirectory("/testem/tracking/");
  fTrackdir->SetGuidance("step length");
  fMaxStepLength = 
        new G4UIcmdWithADoubleAndUnit("/testem/tracking/setMaxStepLength",this);
  fMaxStepLength->SetGuidance("Set the maximum length of tracking step");
  fMaxStepLength->SetGuidance("when integrating magnetic field line.");
  fMaxStepLength->SetParameterName("Size",false);
  fMaxStepLength->SetRange("Size>0.");
  fMaxStepLength->SetUnitCategory("Length");
  fMaxStepLength->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fMaterCmd;
  delete fSizeCmd;
  delete fMaxStepCmd;
  delete fMaxStepLength;
  delete fDetDir;
  delete fTestemDir;
  delete fTrackdir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == fMaterCmd )
   { fDetector->SetMaterial(newValue);}

  if( command == fSizeCmd )
   { fDetector->SetSize(fSizeCmd->GetNewDoubleValue(newValue));}

  if( command == fMaxStepCmd )
   { fDetector->SetMaxStepSize(fMaxStepCmd->GetNewDoubleValue(newValue));}

  if( command == fMaxStepLength )
   { fDetector->SetMaxStepLength(fMaxStepLength->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
