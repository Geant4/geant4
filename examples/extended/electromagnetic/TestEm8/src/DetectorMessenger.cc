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
// $Id: DetectorMessenger.cc,v 1.2 2010-10-27 14:52:07 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/////////////////////////////////////////////////////////////////////////
//
// TestEm8: Gaseous detector
//
// Created: 31.08.2010 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
// 

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

///////////////////////////////////////////////////////////////////////////////////

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:Detector(Det)
{ 
  detDir = new G4UIdirectory("/testem/");
  detDir->SetGuidance("Detector control.");
      
  GasMaterCmd = new G4UIcmdWithAString("/testem/setGasMat",this);
  GasMaterCmd->SetGuidance("Select material of the detector.");
  GasMaterCmd->SetParameterName("gmat",true);
  GasMaterCmd->SetDefaultValue("Argon");
  GasMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  WindowMaterCmd = new G4UIcmdWithAString("/testem/setWindowMat",this);
  WindowMaterCmd->SetGuidance("Select material of the window.");
  WindowMaterCmd->SetParameterName("wmat",true);
  WindowMaterCmd->SetDefaultValue("Mylar");
  WindowMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  WorldMaterCmd = new G4UIcmdWithAString("/testem/setWorldMat",this);
  WorldMaterCmd->SetGuidance("Select material of the world.");
  WorldMaterCmd->SetParameterName("worldmat",true);
  WorldMaterCmd->SetDefaultValue("empty");
  WorldMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  GasThickCmd = new G4UIcmdWithADoubleAndUnit("/testem/setGasThick",this);
  GasThickCmd->SetGuidance("Set thickness of the detector");
  GasThickCmd->SetParameterName("SizeZ",false,false);
  GasThickCmd->SetUnitCategory("Length");
  GasThickCmd->SetDefaultUnit("mm");
  GasThickCmd->SetRange("SizeZ>0.");
  GasThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  GasRadCmd = new G4UIcmdWithADoubleAndUnit("/testem/setGasRad",this);
  GasRadCmd->SetGuidance("Set radius of the detector");
  GasRadCmd->SetParameterName("SizeR",false,false);
  GasRadCmd->SetUnitCategory("Length");
  GasRadCmd->SetDefaultUnit("mm");
  GasRadCmd->SetRange("SizeR>0.");
  GasRadCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  WinThickCmd = new G4UIcmdWithADoubleAndUnit("/testem/setWindowThick",this);
  WinThickCmd->SetGuidance("Set thickness of the window");
  WinThickCmd->SetParameterName("delta",false,false);
  WinThickCmd->SetUnitCategory("Length");
  WinThickCmd->SetDefaultUnit("mm");
  WinThickCmd->SetRange("delta>0.");
  WinThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ionCmd = new G4UIcmdWithADoubleAndUnit("/testem/setPairEnergy",this);
  ionCmd->SetGuidance("Set energy per electron-ion pair for detector");
  ionCmd->SetParameterName("en",false,false);
  ionCmd->SetUnitCategory("Energy");
  ionCmd->SetDefaultUnit("eV");
  ionCmd->SetRange("en>0.");
  ionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
}

//////////////////////////////////////////////////////////////////////////////

DetectorMessenger::~DetectorMessenger()
{
  delete GasMaterCmd; 
  delete GasThickCmd; 
  delete GasRadCmd;  
  delete WinThickCmd; 
  delete WindowMaterCmd;
  delete WorldMaterCmd;
  delete ionCmd; 
  delete detDir;
}

//////////////////////////////////////////////////////////////////////////////////

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == GasMaterCmd )
    { 
      Detector->SetGasMaterial(newValue);
    } 
  else if( command == WindowMaterCmd )
    { 
      Detector->SetContainerMaterial(newValue);
    } 
  else if( command == WorldMaterCmd )
    { 
      Detector->SetWorldMaterial(newValue);
    } 
  else if( command == GasThickCmd )
    { 
      Detector->SetGasThickness(GasThickCmd->GetNewDoubleValue(newValue));
    } 
  else if( command == GasRadCmd )
    { 
      Detector->SetGasRadius(GasRadCmd->GetNewDoubleValue(newValue));
    } 
  else if( command == WinThickCmd )
    { 
      Detector->SetContainerThickness(WinThickCmd->GetNewDoubleValue(newValue));
    }
  else if( command == ionCmd )
    { 
      Detector->SetPairEnergy(ionCmd->GetNewDoubleValue(newValue));
    }
}

/////////////////////////////////////////////////////////////////////////////////
