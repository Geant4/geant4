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
/// \file electromagnetic/MicroElec-SEYv1/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
// $Id: DetectorMessenger.cc 84208 2014-10-10 14:44:50Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithABool.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(),fDetector(Det),
 fTestemDir(0),
 fDetDir(0),    
 fMaterCmd(0),
 fSizeCmd(0)
{ 
  fTestemDir = new G4UIdirectory("/testem/");
  fTestemDir->SetGuidance("commands specific to this example");
  
  fDetDir = new G4UIdirectory("/testem/det/");
  fDetDir->SetGuidance("detector construction");
  
  fMaterCmd = new G4UIcmdWithAString("/testem/det/setMat",this);
  fMaterCmd->SetGuidance("Select material of the box.");
  fMaterCmd->SetParameterName("choice",false);
  fMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fMaterCmd->SetToBeBroadcasted(false);
  
  fMaterSurfCmd = new G4UIcmdWithAString("/testem/det/setMatSurf", this);
  fMaterSurfCmd->SetGuidance("Select material of the box.");
  fMaterSurfCmd->SetParameterName("choice", false);
  fMaterSurfCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fMaterSurfCmd->SetToBeBroadcasted(false);

  fMaterLayer1Cmd = new G4UIcmdWithAString("/testem/det/setMatLayer1", this);
  fMaterLayer1Cmd->SetGuidance("Select material of the box.");
  fMaterLayer1Cmd->SetParameterName("choice", false);
  fMaterLayer1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fMaterLayer1Cmd->SetToBeBroadcasted(false);

  fMaterLayer2Cmd = new G4UIcmdWithAString("/testem/det/setMatLayer2", this);
  fMaterLayer2Cmd->SetGuidance("Select material of the box.");
  fMaterLayer2Cmd->SetParameterName("choice", false);
  fMaterLayer2Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fMaterLayer2Cmd->SetToBeBroadcasted(false);

  fMaterLayer3Cmd = new G4UIcmdWithAString("/testem/det/setMatLayer3", this);
  fMaterLayer3Cmd->SetGuidance("Select material of the box.");
  fMaterLayer3Cmd->SetParameterName("choice", false);
  fMaterLayer3Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fMaterLayer3Cmd->SetToBeBroadcasted(false);

  fMaterLayer4Cmd = new G4UIcmdWithAString("/testem/det/setMatLayer4", this);
  fMaterLayer4Cmd->SetGuidance("Select material of the box.");
  fMaterLayer4Cmd->SetParameterName("choice", false);
  fMaterLayer4Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fMaterLayer4Cmd->SetToBeBroadcasted(false);

  fSizeCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSize",this);
  fSizeCmd->SetGuidance("Set size of the box");
  fSizeCmd->SetParameterName("Size",false);
  fSizeCmd->SetRange("Size>0.");
  fSizeCmd->SetUnitCategory("Length");
  fSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fSizeCmd->SetToBeBroadcasted(false);

  fWidthCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setWidth", this);
  fWidthCmd->SetGuidance("Set Width of the box");
  fWidthCmd->SetParameterName("Width", false);
  fWidthCmd->SetRange("Width>0.");
  fWidthCmd->SetUnitCategory("Length");
  fWidthCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fWidthCmd->SetToBeBroadcasted(false);


  fSizeSurfCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeSurf", this);
  fSizeSurfCmd->SetGuidance("Set size of the box");
  fSizeSurfCmd->SetParameterName("Size", false);
  fSizeSurfCmd->SetRange("Size>0.");
  fSizeSurfCmd->SetUnitCategory("Length");
  fSizeSurfCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSizeSurfCmd->SetToBeBroadcasted(false);

  fSizeLayer1Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeLayer1", this);
  fSizeLayer1Cmd->SetGuidance("Set size of the box");
  fSizeLayer1Cmd->SetParameterName("Size", false);
  fSizeLayer1Cmd->SetRange("Size>0.");
  fSizeLayer1Cmd->SetUnitCategory("Length");
  fSizeLayer1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSizeLayer1Cmd->SetToBeBroadcasted(false);

  fSizeLayer2Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeLayer2", this);
  fSizeLayer2Cmd->SetGuidance("Set size of the box");
  fSizeLayer2Cmd->SetParameterName("Size", false);
  fSizeLayer2Cmd->SetRange("Size>0.");
  fSizeLayer2Cmd->SetUnitCategory("Length");
  fSizeLayer2Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSizeLayer2Cmd->SetToBeBroadcasted(false);

  fSizeLayer3Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeLayer3", this);
  fSizeLayer3Cmd->SetGuidance("Set size of the box");
  fSizeLayer3Cmd->SetParameterName("Size", false);
  fSizeLayer3Cmd->SetRange("Size>0.");
  fSizeLayer3Cmd->SetUnitCategory("Length");
  fSizeLayer3Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSizeLayer3Cmd->SetToBeBroadcasted(false);

  fSizeLayer4Cmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSizeLayer4", this);
  fSizeLayer4Cmd->SetGuidance("Set size of the box");
  fSizeLayer4Cmd->SetParameterName("Size", false);
  fSizeLayer4Cmd->SetRange("Size>0.");
  fSizeLayer4Cmd->SetUnitCategory("Length");
  fSizeLayer4Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSizeLayer4Cmd->SetToBeBroadcasted(false);


  UpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update", this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_PreInit, G4State_Idle);


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fMaterCmd;
  delete fSizeCmd;
  delete fWidthCmd;
  delete fDetDir;  
  delete fTestemDir;
  delete  UpdateCmd;
 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == fMaterCmd )
   { fDetector->SetMaterial(newValue);}

  if (command == fMaterSurfCmd)
  {
	  fDetector->SetMaterialSurface(newValue);
  }
  if (command == fMaterLayer1Cmd)
  {
	  fDetector->SetMaterialLayer1(newValue);
  }
  if (command == fMaterLayer2Cmd)
  {
	  fDetector->SetMaterialLayer2(newValue);
  }
  if (command == fMaterLayer3Cmd)
  {
	  fDetector->SetMaterialLayer3(newValue);
  }
  if (command == fMaterLayer4Cmd)
  {
	  fDetector->SetMaterialLayer4(newValue);
  }


  if( command == fSizeCmd )
   { 
	  fDetector->SetSize(fSizeCmd->GetNewDoubleValue(newValue));
  }
  if (command == fWidthCmd)
  {
	  fDetector->SetWidth(fWidthCmd->GetNewDoubleValue(newValue));
  }
  if (command == fSizeSurfCmd)
  {
	  fDetector->SetSizeSurface(fSizeSurfCmd->GetNewDoubleValue(newValue));
  }
  if (command == fSizeLayer1Cmd)
  {
	  fDetector->SetSizeLayer1(fSizeLayer1Cmd->GetNewDoubleValue(newValue));
  }
  if (command == fSizeLayer2Cmd)
  {
	  fDetector->SetSizeLayer2(fSizeLayer2Cmd->GetNewDoubleValue(newValue));
  }
  if (command == fSizeLayer3Cmd)
  {
	  fDetector->SetSizeLayer3(fSizeLayer3Cmd->GetNewDoubleValue(newValue));
  }
  if (command == fSizeLayer4Cmd)
  {
	  fDetector->SetSizeLayer4(fSizeLayer4Cmd->GetNewDoubleValue(newValue));
  }

  if (command == UpdateCmd) {
	  fDetector->UpdateGeometry();
  }
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
