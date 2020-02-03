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
/// \file electromagnetic/TestEm9/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
//
//
/////////////////////////////////////////////////////////////////////////
//
// TestEm9: Crystal calorimeter
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * det)
  :G4UImessenger(),fDetector(det),
   fAtestemDir(0),
   fAMaterCmd(0),
   fALBinCmd(0),
   fAl1Cmd(0),
   fAl2Cmd(0),
   fAl3Cmd(0),
   fAl4Cmd(0),
   fAl5Cmd(0),
   fAl6Cmd(0),
   fAUpdateCmd(0),
   fAaccCmd1(0),
   fAaccCmd2(0),
   fAaccCmd3(0)
{
  fAtestemDir = new G4UIdirectory("/testem/");
  fAtestemDir->SetGuidance(" detector control.");

  fAMaterCmd = new G4UIcmdWithAString("/testem/det/CalMat",this);
  fAMaterCmd->SetGuidance("Select Material for calorimeter");
  fAMaterCmd->SetParameterName("calMaterial",false);
  fAMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fALBinCmd = new G4UIcmdWithAString("/testem/det/AbsMat",this);
  fALBinCmd->SetGuidance("Select Material for absorber");
  fALBinCmd->SetParameterName("absMarerial",false);
  fALBinCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAl1Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/EcalLength",this);
  fAl1Cmd->SetGuidance("Set length of Ecal");
  fAl1Cmd->SetParameterName("lEcal",false);
  fAl1Cmd->SetUnitCategory("Length");
  fAl1Cmd->SetRange("lEcal>0");
  fAl1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAl2Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/EcalWidth",this);
  fAl2Cmd->SetGuidance("Set width of Ecal crystal");
  fAl2Cmd->SetParameterName("wEcal",false);
  fAl2Cmd->SetUnitCategory("Length");
  fAl2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAl3Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/AbsLength",this);
  fAl3Cmd->SetGuidance("Set length of the absorber");
  fAl3Cmd->SetParameterName("lAbs",false);
  fAl3Cmd->SetUnitCategory("Length");
  fAl3Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAl4Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/VertexLength",this);
  fAl4Cmd->SetGuidance("Set length of the vertex region");
  fAl4Cmd->SetParameterName("lVert",false);
  fAl4Cmd->SetUnitCategory("Length");
  fAl4Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAl5Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/PadLength",this);
  fAl5Cmd->SetGuidance("Set length of vertex detector");
  fAl5Cmd->SetParameterName("lPad",false);
  fAl5Cmd->SetUnitCategory("Length");
  fAl5Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAl6Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/PadWidth",this);
  fAl6Cmd->SetGuidance("Set width of a vertex pad");
  fAl6Cmd->SetParameterName("wPad",false);
  fAl6Cmd->SetUnitCategory("Length");
  fAl6Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAUpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update",this);
  fAUpdateCmd->SetGuidance("Update geometry.");
  fAUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fAUpdateCmd->SetGuidance("if you changed geometrical value(s)");
  fAUpdateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAaccCmd1 = new G4UIcmdWith3Vector("/testem/det/acceptance1",this);
  fAaccCmd1->SetGuidance("set Edep and RMS");
  fAaccCmd1->SetGuidance("acceptance values for central cell");
  fAaccCmd1->SetParameterName("edep","rms","limit",true);
  fAaccCmd1->SetRange("edep>0 && edep<1 && rms>0");
  fAaccCmd1->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAaccCmd2 = new G4UIcmdWith3Vector("/testem/det/acceptance9",this);
  fAaccCmd2->SetGuidance("set Edep and RMS");
  fAaccCmd2->SetGuidance("acceptance values for 3x3 matrix");
  fAaccCmd2->SetParameterName("edep","rms","limit",true);
  fAaccCmd2->SetRange("edep>0 && edep<1 && rms>0");
  fAaccCmd2->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAaccCmd3 = new G4UIcmdWith3Vector("/testem/det/acceptance25",this);
  fAaccCmd3->SetGuidance("set Edep and RMS");
  fAaccCmd3->SetGuidance("acceptance values for 5x5 matrix");
  fAaccCmd3->SetParameterName("edep","rms","limit",true);
  fAaccCmd3->SetRange("edep>0 && edep<1 && rms>0");
  fAaccCmd3->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fAMaterCmd;
  delete fALBinCmd;
  delete fAl1Cmd;
  delete fAl2Cmd;
  delete fAl3Cmd;
  delete fAl4Cmd;
  delete fAl5Cmd;
  delete fAl6Cmd;
  delete fAUpdateCmd;
  delete fAtestemDir;
  delete fAaccCmd1;
  delete fAaccCmd2;
  delete fAaccCmd3;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == fAMaterCmd )
   { fDetector->SetEcalMaterial(newValue);}

  if( command == fALBinCmd )
   { fDetector->SetAbsMaterial(newValue);}

  if( command == fAl1Cmd )
   { fDetector->SetEcalLength(fAl1Cmd->GetNewDoubleValue(newValue));}

  if( command == fAl2Cmd )
   { fDetector->SetEcalWidth(fAl2Cmd->GetNewDoubleValue(newValue));}

  if( command == fAl3Cmd )
   { fDetector->SetAbsLength(fAl3Cmd->GetNewDoubleValue(newValue));}

  if( command == fAl4Cmd )
   { fDetector->SetVertexLength(fAl4Cmd->GetNewDoubleValue(newValue));}

  if( command == fAl5Cmd )
   { fDetector->SetPadLength(fAl5Cmd->GetNewDoubleValue(newValue));}

  if( command == fAl6Cmd )
   { fDetector->SetPadWidth(fAl6Cmd->GetNewDoubleValue(newValue));}

  if( command == fAUpdateCmd )
   { fDetector->UpdateGeometry();}

  HistoManager* histo = HistoManager::GetPointer();
  if( command == fAaccCmd1 )
   { histo->SetEdepAndRMS(0,fAaccCmd1->GetNew3VectorValue(newValue));}

  if( command == fAaccCmd2 )
   { histo->SetEdepAndRMS(1,fAaccCmd2->GetNew3VectorValue(newValue));}

  if( command == fAaccCmd3 )
   { histo->SetEdepAndRMS(2,fAaccCmd3->GetNew3VectorValue(newValue));}


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

