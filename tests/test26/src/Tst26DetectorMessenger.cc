//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: Tst26DetectorMessenger.cc,v 1.3 2003-02-06 11:53:27 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
/////////////////////////////////////////////////////////////////////////
//
// test26: Cut per region physics
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst26DetectorMessenger.hh"

#include "Tst26DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26DetectorMessenger::Tst26DetectorMessenger(Tst26DetectorConstruction * Det)
:Tst26Detector(Det)
{ 
  testemDir = new G4UIdirectory("/testem/");
  testemDir->SetGuidance("Tst26 detector control.");
      
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26DetectorMessenger::~Tst26DetectorMessenger()
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == MaterCmd )
   { Tst26Detector->SetEcalMaterial(newValue);}
   
  if( command == LBinCmd )
   { Tst26Detector->SetAbsMaterial(newValue);}
   
  if( command == l1Cmd )
   { Tst26Detector->SetEcalLength(l1Cmd->GetNewDoubleValue(newValue));}
      
  if( command == l2Cmd )
   { Tst26Detector->SetEcalWidth(l2Cmd->GetNewDoubleValue(newValue));}
     
  if( command == l3Cmd )
   { Tst26Detector->SetAbsLength(l3Cmd->GetNewDoubleValue(newValue));}
      
  if( command == l4Cmd )
   { Tst26Detector->SetVertexLength(l4Cmd->GetNewDoubleValue(newValue));}
     
  if( command == l5Cmd )
   { Tst26Detector->SetPadLength(l5Cmd->GetNewDoubleValue(newValue));}
      
  if( command == l6Cmd )
   { Tst26Detector->SetPadWidth(l6Cmd->GetNewDoubleValue(newValue));}
     
  if( command == UpdateCmd )
   { Tst26Detector->UpdateGeometry();}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

