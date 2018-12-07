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
/// \file field/field04/src/F04DetectorMessenger.cc
/// \brief Implementation of the F04DetectorMessenger class
//

#include "F04DetectorMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04DetectorMessenger::F04DetectorMessenger(F04DetectorConstruction* detector)
 : fDetector(detector)
{
  fDetDir = new G4UIdirectory("/field04/");
  fDetDir->SetGuidance(" field04 Simulation ");

  fWorldMaterCmd = new G4UIcmdWithAString("/field04/SetWorldMat",this);
  fWorldMaterCmd->SetGuidance("Select Material of the World");
  fWorldMaterCmd->SetParameterName("wchoice",true);
  fWorldMaterCmd->SetDefaultValue("Air");
  fWorldMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fWorldMaterCmd->SetToBeBroadcasted(false);

  fWorldRCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetWorldR",this);
  fWorldRCmd->SetGuidance("Set Radius of the World");
  fWorldRCmd->SetParameterName("WSizeR",false,false);
  fWorldRCmd->SetDefaultUnit("cm");
  fWorldRCmd->SetRange("WSizeR>0.");
  fWorldRCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fWorldRCmd->SetToBeBroadcasted(false);

  fWorldZCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetWorldZ",this);
  fWorldZCmd->SetGuidance("Set Length of the World");
  fWorldZCmd->SetParameterName("WSizeZ",false,false);
  fWorldZCmd->SetDefaultUnit("cm");
  fWorldZCmd->SetRange("WSizeZ>0.");
  fWorldZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fWorldZCmd->SetToBeBroadcasted(false);

  fCaptureRCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetCaptureR",this);
  fCaptureRCmd->SetGuidance("Set Radius of the Capture Magnet");
  fCaptureRCmd->SetParameterName("CSizeR",false,false);
  fCaptureRCmd->SetDefaultUnit("cm");
  fCaptureRCmd->SetRange("CSizeR>0.");
  fCaptureRCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fCaptureRCmd->SetToBeBroadcasted(false);

  fCaptureZCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetCaptureZ",this);
  fCaptureZCmd->SetGuidance("Set Length of the Capture Magnet");
  fCaptureZCmd->SetParameterName("CSizeZ",false,false);
  fCaptureZCmd->SetDefaultUnit("cm");
  fCaptureZCmd->SetRange("CSizeZ>0.");
  fCaptureZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fCaptureZCmd->SetToBeBroadcasted(false);

  fTransferRCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetTransferR",this);
  fTransferRCmd->SetGuidance("Set Radius of the Transfer Magnet");
  fTransferRCmd->SetParameterName("TSizeR",false,false);
  fTransferRCmd->SetDefaultUnit("cm");
  fTransferRCmd->SetRange("TSizeR>0.");
  fTransferRCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fTransferRCmd->SetToBeBroadcasted(false);

  fTransferZCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetTransferZ",this);
  fTransferZCmd->SetGuidance("Set Length of the Transfer Magnet");
  fTransferZCmd->SetParameterName("TSizeZ",false,false);
  fTransferZCmd->SetDefaultUnit("cm");
  fTransferZCmd->SetRange("TSizeZ>0.");
  fTransferZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fTransferZCmd->SetToBeBroadcasted(false);

  fTransferPCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetTransferP",this);
  fTransferPCmd->SetGuidance("Set Z pos of the T-Mgnt from end of C-Mgnt");
  fTransferPCmd->SetParameterName("TSizeP",false,false);
  fTransferPCmd->SetDefaultUnit("cm");
  fTransferPCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fTransferPCmd->SetToBeBroadcasted(false);

  fTgtMaterCmd = new G4UIcmdWithAString("/field04/SetTgtMat",this);
  fTgtMaterCmd->SetGuidance("Select Material of the Target");
  fTgtMaterCmd->SetParameterName("tchoice",true);
  fTgtMaterCmd->SetDefaultValue("Tungsten");
  fTgtMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fTgtMaterCmd->SetToBeBroadcasted(false);

  fTgtRadCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetTgtRad",this);
  fTgtRadCmd->SetGuidance("Set Radius of the Target");
  fTgtRadCmd->SetParameterName("TgtSizeR",false,false);
  fTgtRadCmd->SetDefaultUnit("cm");
  fTgtRadCmd->SetRange("TgtSizeR>0.");
  fTgtRadCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fTgtRadCmd->SetToBeBroadcasted(false);

  fTgtThickCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetTgtThick",this);
  fTgtThickCmd->SetGuidance("Set Thickness of the Target");
  fTgtThickCmd->SetParameterName("TgtSizeZ",false,false);
  fTgtThickCmd->SetDefaultUnit("cm");
  fTgtThickCmd->SetRange("TgtSizeZ>0.");
  fTgtThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fTgtThickCmd->SetToBeBroadcasted(false);

  fTgtPosCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetTgtPos",this);
  fTgtPosCmd->SetGuidance("Set Z pos of the tgt relative to C-Mgnt centre");
  fTgtPosCmd->SetParameterName("TgtSizeP",false,false);
  fTgtPosCmd->SetDefaultUnit("cm");
  fTgtPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fTgtPosCmd->SetToBeBroadcasted(false);

  fTgtAngCmd = new G4UIcmdWithAnInteger("/field04/SetTgtAng",this);
  fTgtAngCmd->
    SetGuidance("Set the angle [in deg] of the Tgt relative to C-Mgnt centre");
  fTgtAngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fTgtAngCmd->SetToBeBroadcasted(false);

  fDgrMaterCmd = new G4UIcmdWithAString("/field04/SetDgrMat",this);
  fDgrMaterCmd->SetGuidance("Select Material of the Degrader");
  fDgrMaterCmd->SetParameterName("dchoice",true);
  fDgrMaterCmd->SetDefaultValue("Lead");
  fDgrMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fDgrMaterCmd->SetToBeBroadcasted(false);

  fDgrRadCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetDgrRad",this);
  fDgrRadCmd->SetGuidance("Set Radius of the Degrader");
  fDgrRadCmd->SetParameterName("DrgSizeR",false,false);
  fDgrRadCmd->SetDefaultUnit("cm");
  fDgrRadCmd->SetRange("DrgSizeR>0.");
  fDgrRadCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fDgrRadCmd->SetToBeBroadcasted(false);

  fDgrThickCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetDgrThick",this);
  fDgrThickCmd->SetGuidance("Set Thickness of the Degrader");
  fDgrThickCmd->SetParameterName("DgrSizeZ",false,false);
  fDgrThickCmd->SetDefaultUnit("cm");
  fDgrThickCmd->SetRange("DgrSizeZ>0.");
  fDgrThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fDgrThickCmd->SetToBeBroadcasted(false);

  fDgrPosCmd = new G4UIcmdWithADoubleAndUnit("/field04/SetDgrPos",this);
  fDgrPosCmd->SetGuidance("Set Z pos of the Dgr relative to T-Mgnt centre");
  fDgrPosCmd->SetParameterName("DgrSizeP",false,false);
  fDgrPosCmd->SetDefaultUnit("cm");
  fDgrPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fDgrPosCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04DetectorMessenger::~F04DetectorMessenger()
{
  delete fDetDir;

  delete fWorldMaterCmd;
  delete fWorldRCmd;
  delete fWorldZCmd;

  delete fCaptureRCmd;
  delete fCaptureZCmd;

  delete fTransferRCmd;
  delete fTransferZCmd;
  delete fTransferPCmd;

  delete fTgtMaterCmd;
  delete fTgtRadCmd;
  delete fTgtThickCmd;
  delete fTgtPosCmd;
  delete fTgtAngCmd;

  delete fDgrMaterCmd;
  delete fDgrRadCmd;
  delete fDgrThickCmd;
  delete fDgrPosCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == fWorldMaterCmd )
   { fDetector->SetWorldMaterial(newValue);}

  if( command == fTgtMaterCmd )
   { fDetector->SetTargetMaterial(newValue);}
 
  if( command == fDgrMaterCmd )
   { fDetector->SetDegraderMaterial(newValue);}

  if( command == fWorldRCmd )
   { fDetector->SetWorldSizeR(fWorldRCmd->GetNewDoubleValue(newValue));}

  if( command == fWorldZCmd )
   { fDetector->SetWorldSizeZ(fWorldZCmd->GetNewDoubleValue(newValue));}

  if( command == fCaptureRCmd )
    fDetector->SetCaptureMgntRadius(fCaptureRCmd->GetNewDoubleValue(newValue));

  if( command == fCaptureZCmd )
    fDetector->SetCaptureMgntLength(fCaptureZCmd->GetNewDoubleValue(newValue));

  if( command == fTransferRCmd )
   fDetector->SetTransferMgntRadius(fTransferRCmd->GetNewDoubleValue(newValue));

  if( command == fTransferZCmd )
   fDetector->SetTransferMgntLength(fTransferZCmd->GetNewDoubleValue(newValue));

  if( command == fTransferPCmd )
    fDetector->SetTransferMgntPos(fTransferPCmd->GetNewDoubleValue(newValue));

  if( command == fTgtRadCmd )
    fDetector->SetTargetRadius(fTgtRadCmd->GetNewDoubleValue(newValue));

  if( command == fTgtThickCmd )
    fDetector->SetTargetThickness(fTgtThickCmd->GetNewDoubleValue(newValue));

  if( command == fTgtPosCmd )
    fDetector->SetTargetPos(fTgtPosCmd->GetNewDoubleValue(newValue));

  if( command == fTgtAngCmd )
    fDetector->SetTargetAngle(fTgtAngCmd->GetNewIntValue(newValue));

  if( command == fDgrRadCmd )
    fDetector->SetDegraderRadius(fDgrRadCmd->GetNewDoubleValue(newValue));
 
  if( command == fDgrThickCmd )
    fDetector->SetDegraderThickness(fDgrThickCmd->GetNewDoubleValue(newValue));
 
  if( command == fDgrPosCmd )
    fDetector->SetDegraderPos(fDgrPosCmd->GetNewDoubleValue(newValue));
 
  if( command == fWorldZCmd )
    fDetector->SetWorldSizeZ(fWorldZCmd->GetNewDoubleValue(newValue));

  if( command == fWorldRCmd )
    fDetector->SetWorldSizeR(fWorldRCmd->GetNewDoubleValue(newValue));
}
