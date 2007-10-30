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
//

#include "F04DetectorMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

F04DetectorMessenger::F04DetectorMessenger(F04DetectorConstruction * Det)
 : Detector(Det)
{ 
  detDir = new G4UIdirectory("/field04/");
  detDir->SetGuidance(" field04 Simulation ");

  WorldMaterCmd = new G4UIcmdWithAString("/field04/SetWorldMat",this);
  WorldMaterCmd->SetGuidance("Select Material of the World");
  WorldMaterCmd->SetParameterName("wchoice",true);
  WorldMaterCmd->SetDefaultValue("Air");
  WorldMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  WorldRCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetWorldR",this);
  WorldRCmd->SetGuidance("Set Radius of the World");
  WorldRCmd->SetParameterName("WSizeR",false,false);
  WorldRCmd->SetDefaultUnit("cm");
  WorldRCmd->SetRange("WSizeR>0.");
  WorldRCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  WorldZCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetWorldZ",this);
  WorldZCmd->SetGuidance("Set Length of the World");
  WorldZCmd->SetParameterName("WSizeZ",false,false);
  WorldZCmd->SetDefaultUnit("cm");
  WorldZCmd->SetRange("WSizeZ>0.");
  WorldZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  CaptureRCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetCaptureR",this);
  CaptureRCmd->SetGuidance("Set Radius of the Capture Magnet");
  CaptureRCmd->SetParameterName("CSizeR",false,false);
  CaptureRCmd->SetDefaultUnit("cm");
  CaptureRCmd->SetRange("CSizeR>0.");
  CaptureRCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  CaptureZCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetCaptureZ",this);
  CaptureZCmd->SetGuidance("Set Length of the Capture Magnet");
  CaptureZCmd->SetParameterName("CSizeZ",false,false);
  CaptureZCmd->SetDefaultUnit("cm");
  CaptureZCmd->SetRange("CSizeZ>0.");
  CaptureZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  CaptureB1Cmd = new G4UIcmdWithADoubleAndUnit("/field04/SetCaptureB1",this);
  CaptureB1Cmd->SetGuidance("Set B1 of the Capture Magnet");
  CaptureB1Cmd->SetParameterName("CSizeB1",false,false);
  CaptureB1Cmd->SetDefaultUnit("tesla");
  CaptureB1Cmd->SetRange("CSizeB1>0.");
  CaptureB1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  CaptureB2Cmd = new G4UIcmdWithADoubleAndUnit("/field04/SetCaptureB2",this);
  CaptureB2Cmd->SetGuidance("Set B2 of the Capture Magnet");
  CaptureB2Cmd->SetParameterName("CSizeB2",false,false);
  CaptureB2Cmd->SetDefaultUnit("tesla");
  CaptureB2Cmd->SetRange("CSizeB2>0.");
  CaptureB2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TransferRCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetTransferR",this);
  TransferRCmd->SetGuidance("Set Radius of the Transfer Magnet");
  TransferRCmd->SetParameterName("TSizeR",false,false);
  TransferRCmd->SetDefaultUnit("cm");
  TransferRCmd->SetRange("TSizeR>0.");
  TransferRCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TransferZCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetTransferZ",this);
  TransferZCmd->SetGuidance("Set Length of the Transfer Magnet");
  TransferZCmd->SetParameterName("TSizeZ",false,false);
  TransferZCmd->SetDefaultUnit("cm");
  TransferZCmd->SetRange("TSizeZ>0.");
  TransferZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TransferBCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetTransferB",this);
  TransferBCmd->SetGuidance("Set B of the Transfer Magnet");
  TransferBCmd->SetParameterName("TSizeB",false,false);
  TransferBCmd->SetDefaultUnit("tesla");
  TransferBCmd->SetRange("TSizeB>0.");
  TransferBCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TransferPCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetTransferP",this);
  TransferPCmd->SetGuidance("Set Z pos of the T-Mgnt from end of C-Mgnt");
  TransferPCmd->SetParameterName("TSizeP",false,false);
  TransferPCmd->SetDefaultUnit("cm");
  TransferPCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TgtMaterCmd = new G4UIcmdWithAString("/field04/SetTgtMat",this);
  TgtMaterCmd->SetGuidance("Select Material of the Target");
  TgtMaterCmd->SetParameterName("tchoice",true);
  TgtMaterCmd->SetDefaultValue("Tungsten");
  TgtMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TgtRadCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetTgtRad",this);
  TgtRadCmd->SetGuidance("Set Radius of the Target");
  TgtRadCmd->SetParameterName("TgtSizeR",false,false);
  TgtRadCmd->SetDefaultUnit("cm");
  TgtRadCmd->SetRange("TgtSizeR>0.");
  TgtRadCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TgtThickCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetTgtThick",this);
  TgtThickCmd->SetGuidance("Set Thickness of the Target");
  TgtThickCmd->SetParameterName("TgtSizeZ",false,false);
  TgtThickCmd->SetDefaultUnit("cm");
  TgtThickCmd->SetRange("TgtSizeZ>0.");
  TgtThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TgtPosCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetTgtPos",this);
  TgtPosCmd->SetGuidance("Set Z pos of the tgt relative to C-Mgnt centre");
  TgtPosCmd->SetParameterName("TgtSizeP",false,false);
  TgtPosCmd->SetDefaultUnit("cm");
  TgtPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TgtAngCmd = new G4UIcmdWithAnInteger("/field04/SetTgtAng",this);
  TgtAngCmd->SetGuidance("Set the angle [in deg] of the Tgt relative to C-Mgnt centre");
  TgtAngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DgrMaterCmd = new G4UIcmdWithAString("/field04/SetDgrMat",this);
  DgrMaterCmd->SetGuidance("Select Material of the Degrader");
  DgrMaterCmd->SetParameterName("dchoice",true);
  DgrMaterCmd->SetDefaultValue("Lead");
  DgrMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DgrRadCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetDgrRad",this);
  DgrRadCmd->SetGuidance("Set Radius of the Degrader");
  DgrRadCmd->SetParameterName("DrgSizeR",false,false);
  DgrRadCmd->SetDefaultUnit("cm");
  DgrRadCmd->SetRange("DrgSizeR>0.");
  DgrRadCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DgrThickCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetDgrThick",this);
  DgrThickCmd->SetGuidance("Set Thickness of the Degrader");
  DgrThickCmd->SetParameterName("DgrSizeZ",false,false);
  DgrThickCmd->SetDefaultUnit("cm");
  DgrThickCmd->SetRange("DgrSizeZ>0.");
  DgrThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  DgrPosCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetDgrPos",this);
  DgrPosCmd->SetGuidance("Set Z pos of the Dgr relative to T-Mgnt centre");
  DgrPosCmd->SetParameterName("DgrSizeP",false,false);
  DgrPosCmd->SetDefaultUnit("cm");
  DgrPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  UpdateCmd = new G4UIcmdWithoutParameter("/field04/Update",this);
  UpdateCmd->SetGuidance("Update field04 geometry");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

F04DetectorMessenger::~F04DetectorMessenger()
{
  delete detDir;

  delete WorldMaterCmd;
  delete WorldRCmd;
  delete WorldZCmd;

  delete CaptureRCmd;
  delete CaptureZCmd;
  delete CaptureB1Cmd;
  delete CaptureB2Cmd;

  delete TransferRCmd;
  delete TransferZCmd;
  delete TransferBCmd;
  delete TransferPCmd;

  delete TgtMaterCmd;
  delete TgtRadCmd;
  delete TgtThickCmd;
  delete TgtPosCmd;
  delete TgtAngCmd;

  delete DgrMaterCmd;
  delete DgrRadCmd; 
  delete DgrThickCmd;
  delete DgrPosCmd;

  delete UpdateCmd;
}

void F04DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == WorldMaterCmd )
   { Detector->SetWorldMaterial(newValue);}

  if( command == TgtMaterCmd )
   { Detector->SetTargetMaterial(newValue);}
 
  if( command == DgrMaterCmd )
   { Detector->SetDegraderMaterial(newValue);}

  if( command == WorldRCmd )
   { Detector->SetWorldSizeR(WorldRCmd->GetNewDoubleValue(newValue));}

  if( command == WorldZCmd )
   { Detector->SetWorldSizeZ(WorldZCmd->GetNewDoubleValue(newValue));}

  if( command == CaptureRCmd )
   { Detector->SetCaptureMgntRadius(CaptureRCmd->GetNewDoubleValue(newValue));}

  if( command == CaptureZCmd )
   { Detector->SetCaptureMgntLength(CaptureZCmd->GetNewDoubleValue(newValue));}

  if( command == CaptureB1Cmd )
   { Detector->SetCaptureMgntB1(CaptureB1Cmd->GetNewDoubleValue(newValue));}

  if( command == CaptureB2Cmd )
   { Detector->SetCaptureMgntB2(CaptureB2Cmd->GetNewDoubleValue(newValue));}

  if( command == TransferRCmd )
  { Detector->SetTransferMgntRadius(TransferRCmd->GetNewDoubleValue(newValue));}

  if( command == TransferZCmd )
  { Detector->SetTransferMgntLength(TransferZCmd->GetNewDoubleValue(newValue));}

  if( command == TransferBCmd )
   { Detector->SetTransferMgntB(TransferBCmd->GetNewDoubleValue(newValue));}

  if( command == TransferPCmd )
  { Detector->SetTransferMgntPos(TransferPCmd->GetNewDoubleValue(newValue));}

  if( command == TgtRadCmd )
   { Detector->SetTargetRadius(TgtRadCmd->GetNewDoubleValue(newValue));}

  if( command == TgtThickCmd )
   { Detector->SetTargetThickness(TgtThickCmd->GetNewDoubleValue(newValue));}

  if( command == TgtPosCmd )
   { Detector->SetTargetPos(TgtPosCmd->GetNewDoubleValue(newValue));}

  if( command == TgtAngCmd )
   { Detector->SetTargetAngle(TgtAngCmd->GetNewIntValue(newValue));}

  if( command == DgrRadCmd )
   { Detector->SetDegraderRadius(DgrRadCmd->GetNewDoubleValue(newValue));}
 
  if( command == DgrThickCmd )
   { Detector->SetDegraderThickness(DgrThickCmd->GetNewDoubleValue(newValue));}
   
  if( command == DgrPosCmd )
   { Detector->SetDegraderPos(DgrPosCmd->GetNewDoubleValue(newValue));}
   
  if( command == WorldZCmd )
   { Detector->SetWorldSizeZ(WorldZCmd->GetNewDoubleValue(newValue));}
   
  if( command == WorldRCmd )
   { Detector->SetWorldSizeR(WorldRCmd->GetNewDoubleValue(newValue));}
   
  if( command == UpdateCmd )
   { Detector->UpdateGeometry(); }

}
