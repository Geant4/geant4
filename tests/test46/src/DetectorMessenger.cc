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
// $Id: DetectorMessenger.cc,v 1.4 2009-04-06 12:44:16 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/////////////////////////////////////////////////////////////////////////
//
// DetectorMessenger
//
// Created: 13.07.08 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:Detector(Det)
{
  testDir = new G4UIdirectory("/testCalorim/");
  testDir->SetGuidance(" Hadronic Extended Example.");

  matCmd = new G4UIcmdWithAString("/testCalorim/ecalMaterial",this);
  matCmd->SetGuidance("Select Material for Ecal");
  matCmd->SetParameterName("EcalMaterial",false);
  matCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  mat1Cmd = new G4UIcmdWithAString("/testCalorim/worldMaterial",this);
  mat1Cmd->SetGuidance("Select Material for world");
  mat1Cmd->SetParameterName("wMaterial",false);
  mat1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  rCmd = new G4UIcmdWithADoubleAndUnit("/testCalorim/crystalWidth",this);
  rCmd->SetGuidance("Set width of the Ecal");
  rCmd->SetParameterName("crystalW",false);
  rCmd->SetUnitCategory("Length");
  rCmd->SetRange("crystalW>0"); 
  rCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  lCmd = new G4UIcmdWithADoubleAndUnit("/testCalorim/crystalLength",this);
  lCmd->SetGuidance("Set crystalLength of Ecal");
  lCmd->SetParameterName("length",false);
  lCmd->SetUnitCategory("Length");
  lCmd->SetRange("length>0");
  lCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  gCmd = new G4UIcmdWithADoubleAndUnit("/testCalorim/GapWidth",this);
  gCmd->SetGuidance("Set GapWidth of Ecal");
  gCmd->SetParameterName("length1",false);
  gCmd->SetUnitCategory("Length");
  gCmd->SetRange("length1>0");
  gCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  rCmd1 = new G4UIcmdWithADoubleAndUnit("/testCalorim/hcalWidth",this);
  rCmd1->SetGuidance("Set width of the Hcal");
  rCmd1->SetParameterName("HCalW",false);
  rCmd1->SetUnitCategory("Length");
  rCmd1->SetRange("HCalW>0");
  rCmd1->AvailableForStates(G4State_PreInit,G4State_Idle);

  fCmd = new G4UIcmdWithADoubleAndUnit("/testCalorim/magField",this);
  fCmd->SetGuidance("Set magnetic field along X axis");
  fCmd->SetParameterName("field",false);
  fCmd->SetUnitCategory("Magnetic flux density");
  fCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  updateCmd = new G4UIcmdWithoutParameter("/testCalorim/updateGeometry",this);
  updateCmd->SetGuidance("Update geometry.");
  updateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  updateCmd->SetGuidance("if you changed geometrical value(s)");
  updateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  verbCmd = new G4UIcmdWithAnInteger("/testCalorim/verbose",this);
  verbCmd->SetGuidance("Set verbose for ");
  verbCmd->SetParameterName("verb",false);
  verbCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  binCmd = new G4UIcmdWithAnInteger("/testCalorim/nbins",this);
  binCmd->SetGuidance("Set Nbins for Histo ");
  binCmd->SetParameterName("NBins",false);
  binCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  eCmd = new G4UIcmdWithADoubleAndUnit("/testCalorim/maxE",this);
  eCmd->SetGuidance("Set max Energy ");
  eCmd->SetParameterName("maxE",false);
  eCmd->SetUnitCategory("Energy");
  eCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  facCmd1 = new G4UIcmdWithADouble("/testCalorim/factEcal",this);
  facCmd1->SetGuidance("Set Ecal factor ");
  facCmd1->SetParameterName("fac1",false);
  facCmd1->AvailableForStates(G4State_PreInit,G4State_Idle);

  facCmd2 = new G4UIcmdWithADouble("/testCalorim/factHcal",this);
  facCmd2->SetGuidance("Set Hcal factor ");
  facCmd2->SetParameterName("fac2",false);
  facCmd2->AvailableForStates(G4State_PreInit,G4State_Idle);

  mCmd = new G4UIcmdWithoutParameter("/testCalorim/addPreShower",this);
  mCmd->SetGuidance("Build PreShower ");
  mCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete matCmd;
  delete mat1Cmd;
  delete rCmd;
  delete rCmd1;
  delete lCmd;
  delete facCmd1;
  delete facCmd2;
  delete eCmd;
  delete gCmd;
  delete fCmd;
  delete binCmd;
  delete updateCmd;
  delete testDir;
  delete verbCmd;
  delete mCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == matCmd )
    Detector->SetEcalMaterial(newValue);
  else if( command == mat1Cmd )
    Detector->SetWorldMaterial(newValue);
  else if( command == rCmd ) 
    Detector->SetCrystalWidth(rCmd->GetNewDoubleValue(newValue));
  else if( command == rCmd1 ) 
    Detector->SetHcalWidth(rCmd1->GetNewDoubleValue(newValue));
  else if( command == lCmd ) 
    Detector->SetEcalLength(lCmd->GetNewDoubleValue(newValue));
  else if( command == facCmd1 ) 
    HistoManager::GetPointer()->SetFactor1(facCmd1->GetNewDoubleValue(newValue));
  else if( command == facCmd2 ) 
    HistoManager::GetPointer()->SetFactor2(facCmd2->GetNewDoubleValue(newValue));
  else if( command == gCmd ) 
    Detector->SetGapWidth(gCmd->GetNewDoubleValue(newValue));
  else if( command == fCmd ) 
    Detector->SetMagField(fCmd->GetNewDoubleValue(newValue));
  else if( command == binCmd ) 
    HistoManager::GetPointer()->SetNbins(binCmd->GetNewIntValue(newValue));
  else if( command == eCmd ) 
    HistoManager::GetPointer()->SetMaxEnergy(eCmd->GetNewDoubleValue(newValue));
  else if( command == verbCmd )
    HistoManager::GetPointer()->SetVerbose(verbCmd->GetNewIntValue(newValue));
  else if( command == updateCmd )
    Detector->UpdateGeometry();
  else if( command == mCmd )
    Detector->SetBuildPreShower(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

