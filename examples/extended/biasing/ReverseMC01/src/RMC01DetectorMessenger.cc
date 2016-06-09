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
// $Id: RMC01DetectorMessenger.cc,v 1.2 2009-11-27 14:43:25 ldesorgh Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//////////////////////////////////////////////////////////////
//      Class Name:	RMC01DetectorMessenger
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
//////////////////////////////////////////////////////////////

#include "RMC01DetectorMessenger.hh"

#include "RMC01DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RMC01DetectorMessenger::RMC01DetectorMessenger(
                                           RMC01DetectorConstruction* GeneralDet)
:theDetector(GeneralDet)
{ 
  
  GeneralDir= new G4UIdirectory("/RMC01/");
  GeneralDir->SetGuidance("Control of the Geant4 Reverse Monte Carlo example1");
  
  

  detDir = new G4UIdirectory("/RMC01/geometry/");
  detDir->SetGuidance("Geometry control");
  
  
  SetSensitiveVolumeHeightCmd = new G4UIcmdWithADoubleAndUnit("/RMC01/geometry/SetSensitiveVolumeHeight",this);
  SetSensitiveVolumeHeightCmd->SetGuidance("Set the height of the sensitive cylinder");
  SetSensitiveVolumeHeightCmd->AvailableForStates(G4State_PreInit);
   
  SetSensitiveVolumeRadiusCmd = new G4UIcmdWithADoubleAndUnit("/RMC01/geometry/SetSensitiveVolumeRadius",this);
  SetSensitiveVolumeRadiusCmd->SetGuidance("Set the radius of the sensitive cylinder");
  SetSensitiveVolumeRadiusCmd->AvailableForStates(G4State_PreInit);
  
  SetShieldingThicknessCmd = new G4UIcmdWithADoubleAndUnit("/RMC01/geometry/SetShieldingThickness",this);
  SetShieldingThicknessCmd->SetGuidance("Set the thickness of the Aluminum Shielding sphere");
  SetShieldingThicknessCmd->AvailableForStates(G4State_PreInit);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RMC01DetectorMessenger::~RMC01DetectorMessenger()
{
  if (GeneralDir ) delete GeneralDir;
  if (detDir) delete detDir;
  if (SetSensitiveVolumeHeightCmd) delete  SetSensitiveVolumeHeightCmd;
  if (SetSensitiveVolumeRadiusCmd) delete  SetSensitiveVolumeRadiusCmd;
  if (SetShieldingThicknessCmd) delete  SetShieldingThicknessCmd;
   
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RMC01DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{  if( command == SetSensitiveVolumeHeightCmd ){ 
  	theDetector->SetSensitiveVolumeHeight(SetSensitiveVolumeHeightCmd->GetNewDoubleValue(newValue));  
  }
  
  else if( command == SetSensitiveVolumeRadiusCmd ){ 
  	theDetector->SetSensitiveVolumeRadius(SetSensitiveVolumeRadiusCmd->GetNewDoubleValue(newValue));  
  }
  
  else if( command == SetShieldingThicknessCmd ){ 
  	theDetector->SetShieldingThickness(SetShieldingThicknessCmd->GetNewDoubleValue(newValue));  
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
