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
// $Id: DetectorMessenger.cc,v 1.3 2006-06-29 17:03:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:Detector(Det)
{
  testemDir = new G4UIdirectory("/testem/");
  testemDir->SetGuidance(" detector control.");

  MaterCmd = new G4UIcmdWithAString("/testem/det/CalMat",this);
  MaterCmd->SetGuidance("Select Material for calorimeter");
  MaterCmd->SetParameterName("calMaterial",false);
  MaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  LBinCmd = new G4UIcmdWithAString("/testem/det/AbsMat",this);
  LBinCmd->SetGuidance("Select Material for absorber");
  LBinCmd->SetParameterName("absMarerial",false);
  LBinCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  l1Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/EcalLength",this);
  l1Cmd->SetGuidance("Set length of Ecal");
  l1Cmd->SetParameterName("lEcal",false);
  l1Cmd->SetUnitCategory("Length");
  l1Cmd->SetRange("lEcal>0");
  l1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  l2Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/EcalWidth",this);
  l2Cmd->SetGuidance("Set width of Ecal crystal");
  l2Cmd->SetParameterName("wEcal",false);
  l2Cmd->SetUnitCategory("Length");
  l2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  l3Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/AbsLength",this);
  l3Cmd->SetGuidance("Set length of the absorber");
  l3Cmd->SetParameterName("lAbs",false);
  l3Cmd->SetUnitCategory("Length");
  l3Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  l4Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/VertexLength",this);
  l4Cmd->SetGuidance("Set length of the vertex region");
  l4Cmd->SetParameterName("lVert",false);
  l4Cmd->SetUnitCategory("Length");
  l4Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  l5Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/PadLength",this);
  l5Cmd->SetGuidance("Set length of vertex detector");
  l5Cmd->SetParameterName("lPad",false);
  l5Cmd->SetUnitCategory("Length");
  l5Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  l6Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/PadWidth",this);
  l6Cmd->SetGuidance("Set width of a vertex pad");
  l6Cmd->SetParameterName("wPad",false);
  l6Cmd->SetUnitCategory("Length");
  l6Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update",this);
  UpdateCmd->SetGuidance("Update geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s)");
  UpdateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


  accCmd1 = new G4UIcmdWith3Vector("/testem/det/acceptance1",this);
  accCmd1->SetGuidance("set Edep and RMS");
  accCmd1->SetGuidance("acceptance values for central cell");
  accCmd1->SetParameterName("edep","rms","limit",true);
  accCmd1->SetRange("edep>0 && edep<1 && rms>0");
  accCmd1->AvailableForStates(G4State_PreInit,G4State_Idle);

  accCmd2 = new G4UIcmdWith3Vector("/testem/det/acceptance9",this);
  accCmd2->SetGuidance("set Edep and RMS");
  accCmd2->SetGuidance("acceptance values for 3x3 matrix");
  accCmd2->SetParameterName("edep","rms","limit",true);
  accCmd2->SetRange("edep>0 && edep<1 && rms>0");
  accCmd2->AvailableForStates(G4State_PreInit,G4State_Idle);

  accCmd3 = new G4UIcmdWith3Vector("/testem/det/acceptance25",this);
  accCmd3->SetGuidance("set Edep and RMS");
  accCmd3->SetGuidance("acceptance values for 5x5 matrix");
  accCmd3->SetParameterName("edep","rms","limit",true);
  accCmd3->SetRange("edep>0 && edep<1 && rms>0");
  accCmd3->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete MaterCmd;
  delete LBinCmd;
  delete l1Cmd;
  delete l2Cmd;
  delete l3Cmd;
  delete l4Cmd;
  delete l5Cmd;
  delete l6Cmd;
  delete UpdateCmd;
  delete testemDir;
  delete accCmd1;
  delete accCmd2;
  delete accCmd3;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == MaterCmd )
   { Detector->SetEcalMaterial(newValue);}

  if( command == LBinCmd )
   { Detector->SetAbsMaterial(newValue);}

  if( command == l1Cmd )
   { Detector->SetEcalLength(l1Cmd->GetNewDoubleValue(newValue));}

  if( command == l2Cmd )
   { Detector->SetEcalWidth(l2Cmd->GetNewDoubleValue(newValue));}

  if( command == l3Cmd )
   { Detector->SetAbsLength(l3Cmd->GetNewDoubleValue(newValue));}

  if( command == l4Cmd )
   { Detector->SetVertexLength(l4Cmd->GetNewDoubleValue(newValue));}

  if( command == l5Cmd )
   { Detector->SetPadLength(l5Cmd->GetNewDoubleValue(newValue));}

  if( command == l6Cmd )
   { Detector->SetPadWidth(l6Cmd->GetNewDoubleValue(newValue));}

  if( command == UpdateCmd )
   { Detector->UpdateGeometry();}

  HistoManager* histo = HistoManager::GetPointer();
  if( command == accCmd1 )
   { histo->SetEdepAndRMS(0,accCmd1->GetNew3VectorValue(newValue));}

  if( command == accCmd2 )
   { histo->SetEdepAndRMS(1,accCmd2->GetNew3VectorValue(newValue));}

  if( command == accCmd3 )
   { histo->SetEdepAndRMS(2,accCmd3->GetNew3VectorValue(newValue));}


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

