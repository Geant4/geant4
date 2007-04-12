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
// $Id: Tst52DetectorMessenger.cc,v 1.1 2007-04-12 12:00:17 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// author: Susanna Guatelli (guatelli@ge.infn.it)
// 
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// 
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "Tst52DetectorMessenger.hh"
#include "Tst52DetectorConstruction.hh"

Tst52DetectorMessenger::Tst52DetectorMessenger(
                                           Tst52DetectorConstruction* tst50Det)
:detector(tst50Det)
{ 
  tst51Dir = new G4UIdirectory("/target/");
  tst51Dir -> SetGuidance("Tst52 target control.");

  targetMaterialCmd = new G4UIcmdWithAString("/target/setTargetMat",this);
  targetMaterialCmd -> SetGuidance("Select Material of the Target.");
  targetMaterialCmd -> SetParameterName("choice",false);
  targetMaterialCmd -> AvailableForStates(G4State_Idle);
     
  targetThicknessCmd = new G4UIcmdWithADoubleAndUnit("/target/setTargetThick",this);
  targetThicknessCmd -> SetGuidance("Set Thickness of the target");
  targetThicknessCmd -> SetParameterName("Size",false);
  targetThicknessCmd -> SetRange("Size>=0.");
  targetThicknessCmd -> SetUnitCategory("Length");
  targetThicknessCmd -> AvailableForStates(G4State_Idle);
  
  targetXDimensionCmd = new G4UIcmdWithADoubleAndUnit("/target/setTargetX",this);
  targetXDimensionCmd -> SetGuidance("Set X dimension of the target");
  targetXDimensionCmd -> SetParameterName("Size",false);
  targetXDimensionCmd -> SetRange("Size>=0.");
  targetXDimensionCmd -> SetUnitCategory("Length");
  targetXDimensionCmd -> AvailableForStates(G4State_Idle);
  
  targetYDimensionCmd = new G4UIcmdWithADoubleAndUnit("/target/setTargetY",this);
  targetYDimensionCmd -> SetGuidance("Set Y dimension of the target");
  targetYDimensionCmd -> SetParameterName("Size",false);
  targetYDimensionCmd -> SetRange("Size>=0.");
  targetYDimensionCmd -> SetUnitCategory("Length");
  targetYDimensionCmd -> AvailableForStates(G4State_Idle);
  
  UseUserLimitCmd = new G4UIcmdWithABool("/target/setUserLimits",this);
  UseUserLimitCmd -> SetGuidance("switch on/off the UserLimits");
  UseUserLimitCmd -> SetParameterName("choice",false,false);
  UseUserLimitCmd -> AvailableForStates(G4State_Idle);
      
  setStepMaxCmd = new G4UIcmdWithADoubleAndUnit("/target/setMaxStep",this);
  setStepMaxCmd -> SetGuidance("Set the particle Max Step in the target ");
  setStepMaxCmd -> SetParameterName("Size",false);
  setStepMaxCmd -> SetRange("Size>=0.");
  setStepMaxCmd -> SetUnitCategory("Length");
  setStepMaxCmd -> AvailableForStates(G4State_Idle);

  updateCmd = new G4UIcmdWithoutParameter("/target/update",this);
  updateCmd -> SetGuidance("Update calorimeter geometry.");
  updateCmd -> SetGuidance("This command MUST be applied before \"beamOn\" ");
  updateCmd -> SetGuidance("if you changed geometrical value(s).");
  updateCmd -> AvailableForStates(G4State_Idle);

  voxelCmd = new G4UIcmdWithAnInteger("/target/voxelNumber",this);
  voxelCmd -> SetGuidance("Change the number of voxels along the Z axis.");
  voxelCmd -> AvailableForStates(G4State_Idle);
 }

Tst52DetectorMessenger::~Tst52DetectorMessenger()
{
  delete updateCmd;
  delete setStepMaxCmd;
  delete UseUserLimitCmd; 
  delete targetYDimensionCmd; 
  delete targetXDimensionCmd;
  delete targetThicknessCmd; 
  delete targetMaterialCmd; 
  delete tst51Dir;
}

void Tst52DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if(command == targetMaterialCmd)
   {detector -> SetTargetMaterial(newValue);}
   
  if(command == targetThicknessCmd)
   {detector -> SetTargetThickness(targetThicknessCmd
                                               ->GetNewDoubleValue(newValue));}   if(command == targetXDimensionCmd)
   {detector -> SetTargetX(targetXDimensionCmd->GetNewDoubleValue(newValue));}
  
  if(command == targetYDimensionCmd)
   {detector -> SetTargetY(targetYDimensionCmd->GetNewDoubleValue(newValue));}
 
  if(command == UseUserLimitCmd)
   {detector -> SetUserLimits(UseUserLimitCmd->GetNewBoolValue(newValue));}

  if(command == setStepMaxCmd)
   {detector -> SetMaxStepInTarget(setStepMaxCmd->GetNewDoubleValue(newValue));}
  
  if(command == updateCmd)
   {detector -> UpdateGeometry();}

 if(command ==  voxelCmd)
   {detector -> SetVoxelNumber(voxelCmd -> GetNewIntValue(newValue));}
}
