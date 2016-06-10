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
/// \file hadronic/Hadr00/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
//
// $Id: DetectorMessenger.cc 77210 2013-11-22 01:58:38Z adotti $
//
//
/////////////////////////////////////////////////////////////////////////
//
// DetectorMessenger
//
// Created: 20.06.08 V.Ivanchenko
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
#include "G4UIcmdWithoutParameter.hh"
#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(), fDetector(Det)
{
  ftestDir = new G4UIdirectory("/testhadr/");
  ftestDir->SetGuidance(" Hadronic Extended Example.");

  fmatCmd = new G4UIcmdWithAString("/testhadr/TargetMat",this);
  fmatCmd->SetGuidance("Select Material for the target");
  fmatCmd->SetParameterName("tMaterial",false);
  fmatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fmatCmd->SetToBeBroadcasted(false);

  fmat1Cmd = new G4UIcmdWithAString("/testhadr/WorldMat",this);
  fmat1Cmd->SetGuidance("Select Material for world");
  fmat1Cmd->SetParameterName("wMaterial",false);
  fmat1Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fmat1Cmd->SetToBeBroadcasted(false);

  frCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/TargetRadius",this);
  frCmd->SetGuidance("Set radius of the target");
  frCmd->SetParameterName("radius",false);
  frCmd->SetUnitCategory("Length");
  frCmd->SetRange("radius>0");
  frCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  frCmd->SetToBeBroadcasted(false);

  flCmd = new G4UIcmdWithADoubleAndUnit("/testhadr/TargetLength",this);
  flCmd->SetGuidance("Set length of the target");
  flCmd->SetParameterName("length",false);
  flCmd->SetUnitCategory("Length");
  flCmd->SetRange("length>0");
  flCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  flCmd->SetToBeBroadcasted(false);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fmatCmd;
  delete fmat1Cmd;
  delete frCmd;
  delete flCmd;
  delete ftestDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if( command == fmatCmd ) {
    fDetector->SetTargetMaterial(newValue);
  } else if( command == fmat1Cmd ) {
    fDetector->SetWorldMaterial(newValue);
  } else if( command == frCmd ) {
    fDetector->SetTargetRadius(frCmd->GetNewDoubleValue(newValue));
  } else if( command == flCmd ) { 
    fDetector->SetTargetLength(flCmd->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

