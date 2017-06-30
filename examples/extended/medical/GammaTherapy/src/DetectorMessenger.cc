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
// $Id: DetectorMessenger.cc 103469 2017-04-11 07:29:36Z gcosmo $
//
/// \file medical/GammaTherapy/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
// -------------------------------------------------------------
//      GEANT4  test  IBREM
//
// Authors: V.Grichine, V.Ivanchenko
//
// Modified:
//
// 18-02-03 V.Ivanchenko create
//
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//#include "g4root.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* h):
  fDetector(h)
{
  fDetDir = new G4UIdirectory("/testem/");
  fDetDir->SetGuidance("General  commands");
  fDetDir1= new G4UIdirectory("/testem/physics/");
  fDetDir1->SetGuidance(" commands to define physics");
  fDetDir2= new G4UIdirectory("/testem/gun/");
  fDetDir2->SetGuidance(" commands to define gun");

  fAbsMaterCmd = new G4UIcmdWithAString("/testem/target1Material",this);
  fAbsMaterCmd->SetGuidance("Select Material of the target1.");
  fAbsMaterCmd->SetParameterName("Material1",false);
  fAbsMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fWorldMaterCmd = new G4UIcmdWithAString("/testem/target2Material",this);
  fWorldMaterCmd->SetGuidance("Select Material of the target2.");
  fWorldMaterCmd->SetParameterName("Material2",false);
  fWorldMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAbsThickCmd = new G4UIcmdWithADoubleAndUnit("/testem/mylarZ",this);
  fAbsThickCmd->SetGuidance("Set mylarZ");
  fAbsThickCmd->SetParameterName("mylarZ",false);
  fAbsThickCmd->SetUnitCategory("Length");
  fAbsThickCmd->AvailableForStates(G4State_PreInit);

  fAbsGapCmd = new G4UIcmdWithADoubleAndUnit("/testem/delta",this);
  fAbsGapCmd->SetGuidance("Set gap between absorbers");
  fAbsGapCmd->SetParameterName("delta",false);
  fAbsGapCmd->SetUnitCategory("Length");
  fAbsGapCmd->AvailableForStates(G4State_PreInit);

  fAbsSizYZCmd = new G4UIcmdWithADoubleAndUnit("/testem/target1Z",this);
  fAbsSizYZCmd->SetGuidance("Set targeet1Z");
  fAbsSizYZCmd->SetParameterName("target1Z",false);
  fAbsSizYZCmd->SetUnitCategory("Length");
  fAbsSizYZCmd->AvailableForStates(G4State_PreInit);

  fWorldXCmd = new G4UIcmdWithADoubleAndUnit("/testem/target2Z",this);
  fWorldXCmd->SetGuidance("Set target2Z");
  fWorldXCmd->SetParameterName("target2Z",false);
  fWorldXCmd->SetUnitCategory("Length");
  fWorldXCmd->AvailableForStates(G4State_PreInit);

  fXMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/testem/checkShiftZ",this);
  fXMagFieldCmd->SetGuidance("Define checkShftZ");
  fXMagFieldCmd->SetParameterName("CheckSZ",false);
  fXMagFieldCmd->SetUnitCategory("Length");
  fXMagFieldCmd->AvailableForStates(G4State_PreInit);

  fYMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/testem/sdZ",this);
  fYMagFieldCmd->SetGuidance("Define sensitive detector Z");
  fYMagFieldCmd->SetParameterName("sdZ",false);
  fYMagFieldCmd->SetUnitCategory("Length");
  fYMagFieldCmd->AvailableForStates(G4State_PreInit);

  fZMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/testem/sdShiftZ",this);
  fZMagFieldCmd->SetGuidance("Define shift of sensitive detector");
  fZMagFieldCmd->SetParameterName("sdShiftZ",false);
  fZMagFieldCmd->SetUnitCategory("Length");
  fZMagFieldCmd->AvailableForStates(G4State_PreInit);

  fNumOfAbsCmd = new G4UIcmdWithAnInteger("/testem/numberDivR",this);
  fNumOfAbsCmd->SetGuidance("Set number divisions R");
  fNumOfAbsCmd->SetParameterName("NR",false);
  fNumOfAbsCmd->AvailableForStates(G4State_PreInit);

  fNumOfEvt = new G4UIcmdWithAnInteger("/testem/numberDivZ",this);
  fNumOfEvt->SetGuidance("Set number of divisions Z");
  fNumOfEvt->SetParameterName("NZ",false);
  fNumOfEvt->AvailableForStates(G4State_PreInit);

  fVerbCmd = new G4UIcmdWithAnInteger("/testem/verbose",this);
  fVerbCmd->SetGuidance("Set verbose for ");
  fVerbCmd->SetParameterName("verb",false);
  fVerbCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fIntCmd = new G4UIcmdWithAnInteger("/testem/numberDivE",this);
  fIntCmd->SetGuidance("Set number of divisions E");
  fIntCmd->SetParameterName("NZ",false);
  fIntCmd->AvailableForStates(G4State_PreInit);

  fDeltaECmd = new G4UIcmdWithADoubleAndUnit("/testem/maxEnergy",this);
  fDeltaECmd->SetGuidance("Define scale of secondary energy histogram");
  fDeltaECmd->SetParameterName("DeltaE",false);
  fDeltaECmd->SetUnitCategory("Energy");
  fDeltaECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fDetDir;
  delete fDetDir1;
  delete fDetDir2;

  delete fAbsMaterCmd;
  delete fAbsThickCmd;
  delete fAbsGapCmd;
  delete fAbsSizYZCmd;
  delete fWorldMaterCmd;
  delete fWorldXCmd;
  delete fXMagFieldCmd;
  delete fYMagFieldCmd;
  delete fZMagFieldCmd;
  delete fNumOfAbsCmd;
  delete fNumOfEvt;
  delete fVerbCmd;
  delete fIntCmd;
  delete fDeltaECmd;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{

  if( command == fAbsMaterCmd )
   { fDetector->SetTarget1Material(newValue);}

  if( command == fWorldMaterCmd )
   { fDetector->SetTarget2Material(newValue);}

  if( command == fAbsThickCmd )
   { fDetector->SetMylarZ(fAbsThickCmd->GetNewDoubleValue(newValue));}

  if( command == fAbsGapCmd )
   { fDetector->SetGap(fAbsGapCmd->GetNewDoubleValue(newValue));}

  if( command == fAbsSizYZCmd )
   { fDetector->SetTarget1Z(fAbsSizYZCmd->GetNewDoubleValue(newValue));}

  if( command == fWorldXCmd )
   { fDetector->SetTarget2Z(fWorldXCmd->GetNewDoubleValue(newValue));}

  if( command == fXMagFieldCmd )
   { fDetector->SetCheckShiftZ(fXMagFieldCmd->GetNewDoubleValue(newValue));}

  if( command == fYMagFieldCmd )
   { G4double x = fYMagFieldCmd->GetNewDoubleValue(newValue);
     fDetector->SetAbsorberZ(x);
   }

  if( command == fZMagFieldCmd )
   { fDetector->SetAbsorberShiftZ(fZMagFieldCmd->GetNewDoubleValue(newValue));}

  if( command == fNumOfAbsCmd )
   {
     fDetector->SetNumberDivR(fNumOfAbsCmd->GetNewIntValue(newValue));
   }

  if( command == fNumOfEvt )
   { fDetector->SetNumberDivZ(fNumOfEvt->GetNewIntValue(newValue));}

  if( command == fIntCmd )
   { fDetector->SetNumberDivE(fIntCmd->GetNewIntValue(newValue));}
  if( command == fDeltaECmd )
   { fDetector->SetMaxEnergy(fDeltaECmd->GetNewDoubleValue(newValue));}

  if( command == fVerbCmd ){
     G4int ver = fVerbCmd->GetNewIntValue(newValue);
     fDetector->SetVerbose(ver);
   }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
