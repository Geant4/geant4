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
//    **************************************
//    *                                    *
//    *       CellDetectorMessenger.cc     *
//    *                                    *
//    **************************************
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//         Barbara Mascialino (Barbara.Mascialino@ge.infn.it)
// 
// History:
// -----------
// 20 September 2006   S. Guatelli, B. Mascialino   1st implementation
//
// 
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "CellDetectorMessenger.hh"
#include "CellDetectorConstruction.hh"

CellDetectorMessenger::CellDetectorMessenger(
                                           CellDetectorConstruction* cellDet)
:detector(cellDet)
{ 
  cellDir = new G4UIdirectory("/target/");
  cellDir -> SetGuidance("Cell target control.");

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
  UseUserLimitCmd -> AvailableForStates(G4State_PreInit);
      
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
 }

CellDetectorMessenger::~CellDetectorMessenger()
{
  delete updateCmd;
  delete setStepMaxCmd;
  delete UseUserLimitCmd; 
  delete targetYDimensionCmd; 
  delete targetXDimensionCmd;
  delete targetThicknessCmd; 
  delete targetMaterialCmd; 
  delete cellDir;
}

void CellDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
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
}
