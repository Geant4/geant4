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
// Code developed by:
//  S.Guatelli
//
//    *********************************
//    *                               *
//    *    RemSimDetectorMessenger.cc *
//    *                               *
//    *********************************
//
//
// $Id: RemSimDetectorMessenger.cc,v 1.6 2004-05-21 11:14:52 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "RemSimDetectorMessenger.hh"
#include "RemSimDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

RemSimDetectorMessenger::RemSimDetectorMessenger( RemSimDetectorConstruction* Det): detector(Det)
{ 
  vehicleDir = new G4UIdirectory("/configuration/");
  vehicleDir -> SetGuidance("geometry control.");
        
  vehicleCmd = new G4UIcmdWithAString("/configuration/choose",this);
  vehicleCmd -> SetGuidance("Assign the geometrical set-up to G4RunManager."); 
  vehicleCmd -> SetParameterName("choice",true);
  vehicleCmd -> SetCandidates("vehicle moon");
  vehicleCmd -> AvailableForStates(G4State_Idle); 

  shieldingCmd =  new G4UIcmdWithAString("/configuration/AddShielding",this); 
  shieldingCmd -> SetGuidance("Add shielding."); 
  shieldingCmd -> SetParameterName("choice",true);
  shieldingCmd -> SetCandidates("On Off");
  shieldingCmd -> AvailableForStates(G4State_PreInit,G4State_Idle); 

  astronautCmd =  new G4UIcmdWithAString("/configuration/AddAstronaut",this); 
  astronautCmd -> SetGuidance("Add Astronaut."); 
  astronautCmd -> SetParameterName("choice",true);
  astronautCmd -> SetCandidates("On Off");
  astronautCmd -> AvailableForStates(G4State_PreInit,G4State_Idle); 

  // Fix the parameters of the shielding: material and thickness
  shieldingDir = new G4UIdirectory("/shielding/");
  shieldingDir -> SetGuidance("shielding control.");
       
  materialCmd = new G4UIcmdWithAString("/shielding/material",this);
  materialCmd -> SetGuidance("Change the material of the shielding."); 
  materialCmd -> SetParameterName("choice",true);
  materialCmd -> AvailableForStates(G4State_Idle); 

  thicknessCmd =  new G4UIcmdWithADoubleAndUnit("/shielding/thickness",this);
  thicknessCmd -> SetGuidance("Set the thickness of the shielding."); 
  thicknessCmd -> SetParameterName("Size",true);
  thicknessCmd -> SetRange("Size>=0.");
  thicknessCmd -> SetUnitCategory("Length");
  thicknessCmd -> AvailableForStates(G4State_Idle); 
}

RemSimDetectorMessenger::~RemSimDetectorMessenger()
{
  delete thicknessCmd;
  delete materialCmd;
  delete shieldingDir;
  delete astronautCmd; 
  delete shieldingCmd; 
  delete vehicleCmd;
  delete vehicleDir; 
}

void RemSimDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if(command == vehicleCmd)
    detector -> ChooseConfiguration(newValue); 

  if(command == shieldingCmd)
    detector -> AddShielding(newValue); 

  if(command == astronautCmd)
    {
    if (newValue == "On")
    detector -> AddAstronaut(); 
    }
  if(command == materialCmd)
    {
      detector -> ChangeShieldingMaterial(newValue);
    }
 
  if(command == thicknessCmd)
    {
      detector -> ChangeShieldingThickness
                  (thicknessCmd -> GetNewDoubleValue(newValue));
    }
}

