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
// $Id: HadrontherapyDetectorMessenger.cc,v 1.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo, Francesco Di Rosa
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------
#include "HadrontherapyDetectorMessenger.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

// -----------------------------------------------------------------------------
HadrontherapyDetectorMessenger::HadrontherapyDetectorMessenger(
							       HadrontherapyDetectorConstruction* HadrontherapyDet)
  :HadrontherapyDetector(HadrontherapyDet)
{ 
  HadronDir = new G4UIdirectory("/modulator/");
  HadronDir->SetGuidance("Command to rotate the modulator wheel");
  
  detDir = new G4UIdirectory("/modulator/angle/");
  detDir->SetGuidance("Modulator angle control");

  ModulatorAngleCmd = new G4UIcmdWithADoubleAndUnit("/modulator/angle/modAngle",this);
  ModulatorAngleCmd->SetGuidance("Set Modulator Angle");
  ModulatorAngleCmd->SetParameterName("Size",false);
  ModulatorAngleCmd->SetRange("Size>=0.");
  ModulatorAngleCmd->SetUnitCategory("Angle");  
  ModulatorAngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

// -------------------------------------------------------------------------------
HadrontherapyDetectorMessenger::~HadrontherapyDetectorMessenger()
{
  
  delete ModulatorAngleCmd;
  delete detDir;
  delete HadronDir;  
}
// --------------------------------------------------------------------------------
void HadrontherapyDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  if( command == ModulatorAngleCmd )
    { HadrontherapyDetector->SetModulatorAngle(ModulatorAngleCmd->GetNewDoubleValue(newValue));}

}
// ------------------------------------------------------------------------------------
