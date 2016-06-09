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
// $Id$
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
  shieldingCmd -> SetGuidance("Add shielding layer in vehicle configuration."); 
  shieldingCmd -> SetParameterName("choice",true);
  shieldingCmd -> SetCandidates("On Off");
  shieldingCmd -> AvailableForStates(G4State_Idle); 

  SPECmd =  new G4UIcmdWithAString("/configuration/AddSPE",this); 
  SPECmd -> SetGuidance("Add SPE shelter in vehicle configuration."); 
  SPECmd -> SetParameterName("choice",true);
  SPECmd -> SetCandidates("On Off");
  SPECmd -> AvailableForStates(G4State_Idle); 

  roofCmd =  new G4UIcmdWithAString("/configuration/AddRoof",this); 
  roofCmd -> SetGuidance("Add the Roof to the moon habitat."); 
  roofCmd -> SetParameterName("choice",true);
  roofCmd -> SetCandidates("On Off");
  roofCmd -> AvailableForStates(G4State_Idle); 


  // Fix the parameters of the shielding: material and thickness
  shieldingDir = new G4UIdirectory("/shielding/");
  shieldingDir -> SetGuidance("shielding control.");
       
  thicknessCmd =  new G4UIcmdWithADoubleAndUnit("/shielding/thickness",this);
  thicknessCmd -> SetGuidance("Set the thickness of the shielding."); 
  thicknessCmd -> SetParameterName("Size",true);
  thicknessCmd -> SetRange("Size>=0.");
  thicknessCmd -> SetUnitCategory("Length");
  thicknessCmd -> AvailableForStates(G4State_Idle); 

  roofDir = new G4UIdirectory("/roof/");
  roofDir -> SetGuidance("Roof control.");
  thicknessRoofCmd =  new G4UIcmdWithADoubleAndUnit("/roof/thickness",this);
  thicknessRoofCmd -> SetGuidance("Set the thickness of the Roof."); 
  thicknessRoofCmd -> SetParameterName("Size",true);
  thicknessRoofCmd -> SetRange("Size>=0.");
  thicknessRoofCmd -> SetUnitCategory("Length");
  thicknessRoofCmd -> AvailableForStates(G4State_Idle); 
}

RemSimDetectorMessenger::~RemSimDetectorMessenger()
{
  delete thicknessRoofCmd;
  delete roofDir;
  delete thicknessCmd;
  delete shieldingDir;
  delete roofCmd;
  delete SPECmd; 
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

  if(command == SPECmd)
    detector -> AddShelterSPE(newValue); 

  if(command == roofCmd)
    detector -> AddHabitatRoof(newValue); 

  if(command == thicknessCmd)
      detector -> ChangeShieldingThickness
                  (thicknessCmd -> GetNewDoubleValue(newValue));
    
  if(command == thicknessRoofCmd)
      detector -> ChangeRoofThickness
                  (thicknessRoofCmd -> GetNewDoubleValue(newValue));
}

