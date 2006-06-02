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
// $Id: DetectorMessenger.cc,v 1.1 2006-06-02 19:00:00 vnivanch Exp $
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
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
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

  MaterCmd = new G4UIcmdWithAString("/testem/det/PhantomMat",this);
  MaterCmd->SetGuidance("Select Material for phantom");
  MaterCmd->SetParameterName("phMaterial",false);
  MaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  Mat1Cmd = new G4UIcmdWithAString("/testem/det/GapMat",this);
  Mat1Cmd->SetGuidance("Select Material for gap");
  Mat1Cmd->SetParameterName("gMaterial",false);
  Mat1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  Mat2Cmd = new G4UIcmdWithAString("/testem/det/WorldMat",this);
  Mat2Cmd->SetGuidance("Select Material for world");
  Mat2Cmd->SetParameterName("wMaterial",false);
  Mat2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  l1Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/Radius",this);
  l1Cmd->SetGuidance("Set radius of phantom");
  l1Cmd->SetParameterName("radius",false);
  l1Cmd->SetUnitCategory("Length");
  l1Cmd->SetRange("radius>0");
  l1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  l2Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/Width",this);
  l2Cmd->SetGuidance("Set width of phantom slice");
  l2Cmd->SetParameterName("wEcal",false);
  l2Cmd->SetUnitCategory("Length");
  l2Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  l3Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/Gap",this);
  l3Cmd->SetGuidance("Set gap of phantom slice");
  l3Cmd->SetParameterName("wEcal",false);
  l3Cmd->SetUnitCategory("Length");
  l3Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  enCmd = new G4UIcmdWithADoubleAndUnit("/testem/eMaxNeutron",this);
  enCmd->SetGuidance("Set max n energy");
  enCmd->SetParameterName("nEmax",false);
  enCmd->SetUnitCategory("Energy");
  enCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  eeCmd = new G4UIcmdWithADoubleAndUnit("/testem/eMaxElectron",this);
  eeCmd->SetGuidance("Set max e- energy");
  eeCmd->SetParameterName("eEmax",false);
  eeCmd->SetUnitCategory("Energy");
  eeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  egCmd = new G4UIcmdWithADoubleAndUnit("/testem/eMaxGamma",this);
  egCmd->SetGuidance("Set max gamma energy");
  egCmd->SetParameterName("gEmax",false);
  egCmd->SetUnitCategory("Energy");
  egCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NEbinsCmd = new G4UIcmdWithAnInteger("/testem/nBinsE",this);
  NEbinsCmd->SetGuidance("Set number of bins for Energy");
  NEbinsCmd->SetParameterName("NEbins",false);
  NEbinsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NXYbinsCmd = new G4UIcmdWithAnInteger("/testem/nBinsXY",this);
  NXYbinsCmd->SetGuidance("Set number of bins for XY energy distributions");
  NXYbinsCmd->SetParameterName("NXYbins",false);
  NXYbinsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NbinCmd = new G4UIcmdWithAnInteger("/testem/nBinBraggPeak",this);
  NbinCmd->SetGuidance("Set number of bin with Bragg peak");
  NbinCmd->SetParameterName("NbinBr",false);
  NbinCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NumOfAbsCmd = new G4UIcmdWithAnInteger("/testem/det/numberDivZ",this);
  NumOfAbsCmd->SetGuidance("Set number of slices");
  NumOfAbsCmd->SetParameterName("NZ",false);
  NumOfAbsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update",this);
  UpdateCmd->SetGuidance("Update geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s)");
  UpdateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ntupCmd = new G4UIcmdWithABool("/testem/ntuple",this);
  ntupCmd->SetGuidance("Set ntuple to fill");
  ntupCmd->SetParameterName("ntuple",false);
  ntupCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  verbCmd = new G4UIcmdWithAnInteger("/testem/verbose",this);
  verbCmd->SetGuidance("Set verbose for ");
  verbCmd->SetParameterName("verb",false);
  verbCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete MaterCmd;
  delete Mat1Cmd;
  delete Mat2Cmd;
  delete l1Cmd;
  delete l2Cmd;
  delete l3Cmd;
  delete NumOfAbsCmd;
  delete UpdateCmd;
  delete testemDir;
  delete ntupCmd;
  delete verbCmd;
  delete enCmd;
  delete eeCmd;
  delete egCmd;
  delete NEbinsCmd;
  delete NXYbinsCmd;
  delete NbinCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  HistoManager* man = HistoManager::GetPointer();
  if( command == MaterCmd )
   { Detector->SetAbsMaterial(newValue);}

  if( command == Mat1Cmd )
   { Detector->SetGapMaterial(newValue);}

  if( command == Mat2Cmd )
   { Detector->SetWorldMaterial(newValue);}

  if( command == l1Cmd ) {
    G4double x = l1Cmd->GetNewDoubleValue(newValue);
    man->SetAbsRadius(x);
  }

  if( command == l2Cmd ) {
    G4double x = l2Cmd->GetNewDoubleValue(newValue);
    man->SetAbsWidth(x);
  }

  if( command == l3Cmd ) {
    G4double x = l3Cmd->GetNewDoubleValue(newValue);
    man->SetGapWidth(x);
  }

  if( command == enCmd ) {
    G4double e = enCmd->GetNewDoubleValue(newValue);
    man->SetMaxNeutronEnergy(e);
  }

  if( command == eeCmd ) {
    G4double e = eeCmd->GetNewDoubleValue(newValue);
    man->SetMaxElectronEnergy(e);
  }

  if( command == egCmd ) {
    G4double e = egCmd->GetNewDoubleValue(newValue);
    man->SetMaxGammaEnergy(e);
  }
  
  if( command == ntupCmd )
   { man->SetNtuple(ntupCmd->GetNewBoolValue(newValue));}

  if( command == NumOfAbsCmd ) {
    G4int n = NumOfAbsCmd->GetNewIntValue(newValue);
    man->SetNumberOfAbs(n);
  }

  if( command == NEbinsCmd ) {
    G4int n = NEbinsCmd->GetNewIntValue(newValue);
    man->SetNumBinsEnergy(n);
  }

  if( command == NXYbinsCmd ) {
    G4int n = NXYbinsCmd->GetNewIntValue(newValue);
    man->SetNumBinsXY(n);
  }

  if( command == NbinCmd ) {
    G4int n = NbinCmd->GetNewIntValue(newValue);
    man->SetNumBinBragg(n);
  }
  
  if( command == verbCmd ){
     G4int ver = verbCmd->GetNewIntValue(newValue);
     man->SetVerbose(ver);
  }

  if( command == UpdateCmd )
   { Detector->UpdateGeometry();}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

