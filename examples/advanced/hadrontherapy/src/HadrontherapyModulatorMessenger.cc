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
// This is the *BASIC* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// Visit the Hadrontherapy web site (http://www.lns.infn.it/link/Hadrontherapy) to request 
// the *COMPLETE* version of this program, together with its documentation;
// Hadrontherapy (both basic and full version) are supported by the Italian INFN
// Institute in the framework of the MC-INFN Group
//


#include "HadrontherapyModulatorMessenger.hh"
#include "HadrontherapyModulator.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"


    HadrontherapyModulatorMessenger::HadrontherapyModulatorMessenger(HadrontherapyModulator* Mod)
:Modulator(Mod)

{
    modulatorDir = new G4UIdirectory("/modulator/");
    modulatorDir -> SetGuidance("Command to change the modulator wheel proprties");

    modulatorMatCmd = new G4UIcmdWithAString("/modulator/RMWMat",this);
    modulatorMatCmd -> SetGuidance("Set material of modulatorWheel");
    modulatorMatCmd -> SetParameterName("Material",false);
    modulatorMatCmd -> AvailableForStates(G4State_Idle);
    
    modulatorExternalFile=new G4UIcmdWithAString("/modulator/ReadData",this);
    modulatorExternalFile -> SetGuidance("set properties of modulator steps via a external file");
    modulatorExternalFile -> SetParameterName("FileName",true,false);
    modulatorExternalFile -> SetDefaultValue ("default");
    modulatorExternalFile -> AvailableForStates(G4State_Idle);

    modulatorPositionCmd = new G4UIcmdWith3VectorAndUnit("/modulator/position",this);
    modulatorPositionCmd -> SetGuidance("Set position of modulato");
    modulatorPositionCmd -> SetParameterName("PositionAlongX", 
						    "PositionAlongY", 
						    "PositionAlongZ",false);
    modulatorPositionCmd -> SetDefaultUnit("mm");  
    modulatorPositionCmd -> SetUnitCandidates("mm cm m");  
    modulatorPositionCmd -> AvailableForStates(G4State_Idle);

 
    modulatorOuterRadiusCmd = new G4UIcmdWithADoubleAndUnit("/modulator/outRadius",this);
    modulatorOuterRadiusCmd -> SetGuidance("Set size of outer radius");
    modulatorOuterRadiusCmd-> SetParameterName("Size",false);
    modulatorOuterRadiusCmd -> SetDefaultUnit("mm");  
    modulatorOuterRadiusCmd-> SetUnitCandidates("mm cm m");  
    modulatorOuterRadiusCmd-> AvailableForStates(G4State_Idle);

    modulatorInnerRadiusCmd = new G4UIcmdWithADoubleAndUnit("/modulator/innerRadius",this);
    modulatorInnerRadiusCmd -> SetGuidance("Set size of inner radius");
    modulatorInnerRadiusCmd-> SetParameterName("Size",false);
    modulatorInnerRadiusCmd -> SetDefaultUnit("mm");  
    modulatorInnerRadiusCmd-> SetUnitCandidates("mm cm m");  
    modulatorInnerRadiusCmd-> AvailableForStates(G4State_Idle);

    modulatorAngleCmd = new G4UIcmdWithADoubleAndUnit("/modulator/angle",this);
    modulatorAngleCmd -> SetGuidance("Set Modulator Angle");
    modulatorAngleCmd -> SetParameterName("Size",false);
    modulatorAngleCmd -> SetRange("Size>=0.");
    modulatorAngleCmd -> SetUnitCategory("Angle");  
    modulatorAngleCmd -> AvailableForStates(G4State_Idle);
}

 HadrontherapyModulatorMessenger::~ HadrontherapyModulatorMessenger()
{ 
    delete modulatorAngleCmd;  
    delete modulatorMatCmd;  
    delete modulatorPositionCmd; 
    delete modulatorInnerRadiusCmd; 
    delete modulatorOuterRadiusCmd;
    delete modulatorDir;   
}




void  HadrontherapyModulatorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
    if( command == modulatorAngleCmd )
    {Modulator -> SetModulatorAngle
	 (modulatorAngleCmd -> GetNewDoubleValue(newValue));}

   else if( command == modulatorMatCmd )
   {Modulator -> SetModulatorMaterial(newValue);}
   
    else if (command== modulatorExternalFile)
	 {Modulator->GetDataFromFile(newValue);}

   else if( command == modulatorPositionCmd )
    { G4ThreeVector size = modulatorPositionCmd-> GetNew3VectorValue(newValue);
	Modulator -> SetModulatorPosition(size);}

   else if( command == modulatorOuterRadiusCmd )
   { Modulator -> SetModulatorOuterRadius(
	modulatorOuterRadiusCmd -> GetNewDoubleValue(newValue));}

    else if( command == modulatorInnerRadiusCmd )
    { Modulator -> SetModulatorInnerRadius(
	 modulatorInnerRadiusCmd -> GetNewDoubleValue(newValue));}
}

