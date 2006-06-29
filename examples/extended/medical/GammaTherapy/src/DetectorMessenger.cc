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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "Histo.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorMessenger::DetectorMessenger(DetectorConstruction* h):
  hDet(h)
{
  detDir = new G4UIdirectory("/testem/");
  detDir->SetGuidance("General  commands");
  detDir1= new G4UIdirectory("/testem/physics/");
  detDir1->SetGuidance(" commands to define physics");
  detDir2= new G4UIdirectory("/testem/gun/");
  detDir2->SetGuidance(" commands to define gun");

  AbsMaterCmd = new G4UIcmdWithAString("/testem/target1Material",this);
  AbsMaterCmd->SetGuidance("Select Material of the target1.");
  AbsMaterCmd->SetParameterName("Material1",false);
  AbsMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  WorldMaterCmd = new G4UIcmdWithAString("/testem/target2Material",this);
  WorldMaterCmd->SetGuidance("Select Material of the target2.");
  WorldMaterCmd->SetParameterName("Material2",false);
  WorldMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AbsThickCmd = new G4UIcmdWithADoubleAndUnit("/testem/mylarZ",this);
  AbsThickCmd->SetGuidance("Set mylarZ");
  AbsThickCmd->SetParameterName("mylarZ",false);
  AbsThickCmd->SetUnitCategory("Length");
  AbsThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AbsGapCmd = new G4UIcmdWithADoubleAndUnit("/testem/delta",this);
  AbsGapCmd->SetGuidance("Set gap between absorbers");
  AbsGapCmd->SetParameterName("delta",false);
  AbsGapCmd->SetUnitCategory("Length");
  AbsGapCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  AbsSizYZCmd = new G4UIcmdWithADoubleAndUnit("/testem/target1Z",this);
  AbsSizYZCmd->SetGuidance("Set targeet1Z");
  AbsSizYZCmd->SetParameterName("target1Z",false);
  AbsSizYZCmd->SetUnitCategory("Length");
  AbsSizYZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  WorldXCmd = new G4UIcmdWithADoubleAndUnit("/testem/target2Z",this);
  WorldXCmd->SetGuidance("Set target2Z");
  WorldXCmd->SetParameterName("target2Z",false);
  WorldXCmd->SetUnitCategory("Length");
  WorldXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/testem/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  XMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/testem/checkShiftZ",this);
  XMagFieldCmd->SetGuidance("Define checkShftZ");
  XMagFieldCmd->SetParameterName("CheckSZ",false);
  XMagFieldCmd->SetUnitCategory("Length");
  XMagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  YMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/testem/sdZ",this);
  YMagFieldCmd->SetGuidance("Define sensitive detector Z");
  YMagFieldCmd->SetParameterName("sdZ",false);
  YMagFieldCmd->SetUnitCategory("Length");
  YMagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ZMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/testem/sdShiftZ",this);
  ZMagFieldCmd->SetGuidance("Define shift of sensitive detector");
  ZMagFieldCmd->SetParameterName("sdShiftZ",false);
  ZMagFieldCmd->SetUnitCategory("Length");
  ZMagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  HistoCmd = new G4UIcmdWithAString("/testem/histoName",this);
  HistoCmd->SetGuidance("Set the name of the histo file");
  HistoCmd->SetParameterName("histo",false);
  HistoCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  HistoTCmd = new G4UIcmdWithAString("/testem/histoType",this);
  HistoTCmd->SetGuidance("Set the type of the histo file");
  HistoTCmd->SetParameterName("histoT",false);
  HistoTCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ntupCmd = new G4UIcmdWithABool("/testem/ntuple",this);
  ntupCmd->SetGuidance("Set ntuple to fill");
  ntupCmd->SetParameterName("ntuple",false);
  ntupCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NumOfAbsCmd = new G4UIcmdWithAnInteger("/testem/numberDivR",this);
  NumOfAbsCmd->SetGuidance("Set number divisions R");
  NumOfAbsCmd->SetParameterName("NR",false);
  NumOfAbsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NumOfEvt = new G4UIcmdWithAnInteger("/testem/numberDivZ",this);
  NumOfEvt->SetGuidance("Set number of divisions Z");
  NumOfEvt->SetParameterName("NZ",false);
  NumOfEvt->AvailableForStates(G4State_PreInit,G4State_Idle);

  verbCmd = new G4UIcmdWithAnInteger("/testem/verbose",this);
  verbCmd->SetGuidance("Set verbose for ");
  verbCmd->SetParameterName("verb",false);
  verbCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  intCmd = new G4UIcmdWithAnInteger("/testem/numberDivE",this);
  intCmd->SetGuidance("Set number of divisions E");
  intCmd->SetParameterName("NZ",false);
  intCmd->AvailableForStates(G4State_PreInit);

  nhistCmd = new G4UIcmdWithAnInteger("/testem/histoNumber",this);
  nhistCmd->SetGuidance("Set number of histograms to fill");
  nhistCmd->SetParameterName("HistoNumber",false);
  nhistCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  nDebugSCmd = new G4UIcmdWithAnInteger("/testem/nFirstEventToDebug",this);
  nDebugSCmd->SetGuidance("Set number of the first event to debug");
  nDebugSCmd->SetParameterName("nFirstEventToDebug",false);
  nDebugSCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  nDebugECmd = new G4UIcmdWithAnInteger("/testem/nLastEventToDebug",this);
  nDebugECmd->SetGuidance("Set number of the last event to debug");
  nDebugECmd->SetParameterName("nLastEventToDebug",false);
  nDebugECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DeltaECmd = new G4UIcmdWithADoubleAndUnit("/testem/maxEnergy",this);
  DeltaECmd->SetGuidance("Define scale of secondary energy histogram");
  DeltaECmd->SetParameterName("DeltaE",false);
  DeltaECmd->SetUnitCategory("Energy");
  DeltaECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorMessenger::~DetectorMessenger()
{
  delete NumOfAbsCmd;
  delete AbsMaterCmd;
  delete AbsThickCmd;
  delete AbsGapCmd;
  delete AbsSizYZCmd;
  delete WorldMaterCmd;
  delete WorldXCmd;
  delete UpdateCmd;
  delete XMagFieldCmd;
  delete YMagFieldCmd;
  delete ZMagFieldCmd;
  delete HistoCmd;
  delete HistoTCmd;
  delete NumOfEvt;
  delete verbCmd;
  delete intCmd;
  delete nhistCmd;
  delete nDebugSCmd;
  delete nDebugECmd;
  delete detDir;
  delete detDir1;
  delete detDir2;
  delete DeltaECmd;
  delete ntupCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if( command == AbsMaterCmd )
   { hDet->setTarget1Material(newValue);}

  if( command == WorldMaterCmd )
   { hDet->setTarget2Material(newValue);}

  if( command == AbsThickCmd )
   { hDet->setMylarZ(AbsThickCmd->GetNewDoubleValue(newValue));}

  if( command == AbsGapCmd )
   { hDet->setGap(AbsGapCmd->GetNewDoubleValue(newValue));}

  if( command == AbsSizYZCmd )
   { hDet->setTarget1Z(AbsSizYZCmd->GetNewDoubleValue(newValue));}

  if( command == WorldXCmd )
   { hDet->setTarget2Z(WorldXCmd->GetNewDoubleValue(newValue));}

  if( command == UpdateCmd )
   { hDet->UpdateGeometry(); }

  if( command == XMagFieldCmd )
   { hDet->setCheckShiftZ(XMagFieldCmd->GetNewDoubleValue(newValue));}

  if( command == YMagFieldCmd )
   { G4double x = YMagFieldCmd->GetNewDoubleValue(newValue);
     hDet->setAbsorberZ(x);
   }

  if( command == ZMagFieldCmd )
   { hDet->setAbsorberShiftZ(ZMagFieldCmd->GetNewDoubleValue(newValue));}

  if( command == HistoCmd )
   { (Histo::GetPointer())->SetHistoName(newValue);}

  if( command == HistoTCmd )
   { (Histo::GetPointer())->SetHistoType(newValue);}

  if( command == ntupCmd )
   { (Histo::GetPointer())->SetNtuple(ntupCmd->GetNewBoolValue(newValue));}

  if( command == NumOfAbsCmd )
   {(Histo::GetPointer())->SetNumberDivR(NumOfAbsCmd->GetNewIntValue(newValue));}

  if( command == NumOfEvt )
   {(Histo::GetPointer())->SetNumberDivZ(NumOfEvt->GetNewIntValue(newValue));}

  if( command == verbCmd ){
     G4int ver = verbCmd->GetNewIntValue(newValue);
     (Histo::GetPointer())->SetVerbose(ver);
   }

  if( command == intCmd )
   {(Histo::GetPointer())->SetNumberDivE(intCmd->GetNewIntValue(newValue));}

  if( command == nhistCmd )
   { (Histo::GetPointer())->SetHistoNumber(nhistCmd->GetNewIntValue(newValue));}

  if( command == nDebugSCmd )
   {(Histo::GetPointer())->SetFirstEventToDebug(nDebugSCmd->GetNewIntValue(newValue));}

  if( command == nDebugECmd )
   {(Histo::GetPointer())->SetLastEventToDebug(nDebugECmd->GetNewIntValue(newValue));}

  if( command == DeltaECmd )
   {(Histo::GetPointer()) ->SetMaxEnergy(DeltaECmd->GetNewDoubleValue(newValue));}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
