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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy


#include "TrentoPassiveProtonBeamLine.hh"
#include "TrentoPassiveProtonBeamLineMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4SystemOfUnits.hh"

   TrentoPassiveProtonBeamLineMessenger::TrentoPassiveProtonBeamLineMessenger(TrentoPassiveProtonBeamLine* TrentoLine)
:TrentoPassiveProton(TrentoLine)

{
    changeTheBeamLineDir = new G4UIdirectory("/ChangeBeamLine/");
    changeTheBeamLineDir -> SetGuidance("Command to change the transport beam line");

    changeTheBeamLineNameCmd = new G4UIcmdWithAString("/ChangeBeamLine/beamLineName",this);
    changeTheBeamLineNameCmd -> SetGuidance("Insert the name of the beam line you want simulate");
    changeTheBeamLineNameCmd -> SetParameterName("List",false);
    changeTheBeamLineNameCmd -> AvailableForStates(G4State_PreInit); 

    modulatorDir = new G4UIdirectory("/modulator/");
    modulatorDir -> SetGuidance("Command to rotate the modulator wheel");

    beamLineDir = new G4UIdirectory("/beamLine/");
    beamLineDir -> SetGuidance("set specification of range shifter");  

    ScatteringFoilXSizeCmd = new G4UIcmdWithADoubleAndUnit("/beamLine/ScatteringFoil/thickness",this);
    ScatteringFoilXSizeCmd -> SetGuidance("Set half thickness of scattering foil");
    ScatteringFoilXSizeCmd -> SetParameterName("Size",false);
    ScatteringFoilXSizeCmd -> SetDefaultUnit("mm");  
    ScatteringFoilXSizeCmd -> SetUnitCandidates("mm cm m");  
    ScatteringFoilXSizeCmd -> AvailableForStates(G4State_Idle);

    scatteringFoilMatCmd = new G4UIcmdWithAString("/beamLine/ScatteringFoil/material",this);
    scatteringFoilMatCmd -> SetGuidance("Set material of scatterer");
    scatteringFoilMatCmd -> SetParameterName("choice",false);
    scatteringFoilMatCmd -> AvailableForStates(G4State_Idle);
    
    preCollimatorXSizeCmd = new G4UIcmdWithADoubleAndUnit("/beamLine/preCollimator/XSize",this);
    preCollimatorXSizeCmd -> SetGuidance("Set half x side of pre collimator");
    preCollimatorXSizeCmd -> SetParameterName("Size",false);
    preCollimatorXSizeCmd -> SetDefaultUnit("mm");  
    preCollimatorXSizeCmd -> SetUnitCandidates("mm cm m");  
    preCollimatorXSizeCmd -> AvailableForStates(G4State_Idle);

    preCollimatorXPositionCmd = new G4UIcmdWithADoubleAndUnit("/beamLine/preCollimator/XPosition",this);
    preCollimatorXPositionCmd -> SetGuidance("Set the position along x ");
    preCollimatorXPositionCmd -> SetParameterName("Position",false);
    preCollimatorXPositionCmd -> SetDefaultUnit("mm");  
    preCollimatorXPositionCmd -> SetUnitCandidates("mm cm m");  
    preCollimatorXPositionCmd -> AvailableForStates(G4State_Idle);

    AirTubeYSizeCmd = new G4UIcmdWithADoubleAndUnit("/beamLine/airTube/YSide",this);
    AirTubeYSizeCmd -> SetGuidance("Set of pre collimator the outer radius of the hole in the collimator ");
    AirTubeYSizeCmd -> SetParameterName("Radius",false);
    AirTubeYSizeCmd -> SetDefaultUnit("mm");  
    AirTubeYSizeCmd -> SetUnitCandidates("mm cm m");  
    AirTubeYSizeCmd -> AvailableForStates(G4State_Idle);

    
    AirTubeZSizeCmd = new G4UIcmdWithADoubleAndUnit("/beamLine/FirstCollimator/innerRadius",this);
    AirTubeZSizeCmd -> SetGuidance("Set the inner radius of the collimator ");
    AirTubeZSizeCmd -> SetParameterName("Radius",false);
    AirTubeZSizeCmd -> SetDefaultUnit("mm");  
    AirTubeZSizeCmd -> SetUnitCandidates("mm cm m");  
    AirTubeZSizeCmd -> AvailableForStates(G4State_Idle);

}

TrentoPassiveProtonBeamLineMessenger::~TrentoPassiveProtonBeamLineMessenger()
{ 
    delete AirTubeZSizeCmd;
    delete AirTubeYSizeCmd;
    delete ScatteringFoilXSizeCmd;
    delete scatteringFoilMatCmd;
    delete preCollimatorXSizeCmd;
    delete beamLineDir; 
    delete modulatorDir;   
}




void TrentoPassiveProtonBeamLineMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
    if( command == ScatteringFoilXSizeCmd )
    {TrentoPassiveProton -> SetScatteringFoilXSize
	(ScatteringFoilXSizeCmd -> GetNewDoubleValue(newValue));}

     else if( command == scatteringFoilMatCmd )
    { TrentoPassiveProton -> SetScattererMaterial(newValue);}

    else if( command == preCollimatorXSizeCmd )
    { TrentoPassiveProton -> SetPreCollimatorXSize
	(preCollimatorXSizeCmd -> GetNewDoubleValue(newValue));}

    else if( command == preCollimatorXPositionCmd )
    { TrentoPassiveProton -> SetPreCollimatorXPosition
	(preCollimatorXPositionCmd -> GetNewDoubleValue(newValue));}
   
    else if( command == AirTubeYSizeCmd )
    { TrentoPassiveProton -> SetAirTubeYSize
	(AirTubeYSizeCmd -> GetNewDoubleValue(newValue));}

    else if( command == AirTubeZSizeCmd )
    { TrentoPassiveProton -> SetAirTubeZSize
	(AirTubeZSizeCmd -> GetNewDoubleValue(newValue));}

}

