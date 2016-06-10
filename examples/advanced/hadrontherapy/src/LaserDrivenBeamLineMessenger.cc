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

#include "LaserDrivenBeamLineMessenger.hh"
#include "LaserDrivenBeamLine.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4SystemOfUnits.hh"

/////////////////////////////////////////////////////////////////////////////
LaserDrivenBeamLineMessenger::LaserDrivenBeamLineMessenger(LaserDrivenBeamLine* laserDriven)
  :laserDrivenMessengerPointer(laserDriven)
{
// Messenger directories for the Energy Selector module
  laserDrivenDir = new G4UIdirectory("/LaserDriven/");
  laserDrivenDir -> SetGuidance("The Laser Driven Beam Line module of Hadrontherapy");

  // Messenger directories for the Energy Selector module
  energySelectorDir = new G4UIdirectory("/LaserDriven/EnergySelector/");
  energySelectorDir -> SetGuidance("The Energy selector (ESS) module of Hadrontherapy");

  FcollimatorDir =  new G4UIdirectory("/LaserDriven/EnergySelector/FirstCollimator/");
  FcollimatorDir -> SetGuidance("Define geometrical characteristics of the ESS first collimator");

  ScollimatorDir =  new G4UIdirectory("/LaserDriven/EnergySelector/SecondCollimator/");
  ScollimatorDir -> SetGuidance("Define geometrical characteristics of the ESS second collimator");

  slitDir =  new G4UIdirectory("/LaserDriven/EnergySelector/Slit/");
  slitDir -> SetGuidance("Define geometrical characteristics of the ESS slit");

 // Messenger directories for the Quadrupole module
  quadrupoleDir = new G4UIdirectory("/LaserDriven/Quadrupoles/");
  quadrupoleDir -> SetGuidance("The Quadrupoles module of Hadrontherapy");

  relativePosDir = new G4UIdirectory("/LaserDriven/Quadrupoles/xRelPosition/");
  relativePosDir -> SetGuidance("Define the x relative positions of the quadrupoles");

  // ESS DISABLE
  DisableESSCmd = new G4UIcmdWithoutParameter("/LaserDriven/EnergySelector/Disable", this);	
  DisableESSCmd -> SetGuidance("Disable the Energy Selector.");
  DisableESSCmd -> SetGuidance("This command MUST be applied before \"beamOn\" ");
  DisableESSCmd -> AvailableForStates(G4State_Idle);

  // THE FIRST ESS COLLIMATOR
  //
  // Diameter of the first collimator
  FcollimatorRadiusCmd = new G4UIcmdWithADoubleAndUnit("/LaserDriven/EnergySelector/FirstCollimator/Radius", this);
  FcollimatorRadiusCmd -> SetGuidance("Set the Radius of the first collimator");
  FcollimatorRadiusCmd -> SetParameterName("Size",false);
  FcollimatorRadiusCmd -> SetDefaultUnit("mm"); 
  FcollimatorRadiusCmd -> SetUnitCandidates("mm cm m");  
  FcollimatorRadiusCmd -> AvailableForStates(G4State_Idle);
	
  // Thickness of the first collimator
  FcollimatorThicknessCmd = new G4UIcmdWithADoubleAndUnit("/LaserDriven/EnergySelector/FirstCollimator/Thickness", this);
  FcollimatorThicknessCmd -> SetGuidance("Set the thickness of the first collimator");
  FcollimatorThicknessCmd -> SetParameterName("Size",false);
  FcollimatorThicknessCmd -> SetDefaultUnit("mm"); 
  FcollimatorThicknessCmd -> SetUnitCandidates("mm cm m");  
  FcollimatorThicknessCmd -> AvailableForStates(G4State_Idle);

  // Z Position of the first collimator hole 
  FcollimatorZpositionCmd = new G4UIcmdWithADoubleAndUnit("/LaserDriven/EnergySelector/FirstCollimator/zPosizion", this);
  FcollimatorZpositionCmd -> SetGuidance("Set the Z position of the first collimator hole as respect the internal vacuum chamber center axis");
  FcollimatorZpositionCmd -> SetParameterName("Size",false);
  FcollimatorZpositionCmd -> SetDefaultUnit("mm"); 
  FcollimatorZpositionCmd -> SetUnitCandidates("mm cm m");  
  FcollimatorZpositionCmd -> AvailableForStates(G4State_Idle);
		
  // THE SECOND ESS COLLIMATOR
  //
 // Diameter of the second collimator
  ScollimatorRadiusCmd = new G4UIcmdWithADoubleAndUnit("/LaserDriven/EnergySelector/SecondCollimator/Radius", this);
  ScollimatorRadiusCmd -> SetGuidance("Set the Radius of the second collimator");
  ScollimatorRadiusCmd -> SetParameterName("Size",false);
  ScollimatorRadiusCmd -> SetDefaultUnit("mm"); 
  ScollimatorRadiusCmd -> SetUnitCandidates("mm cm m");  
  ScollimatorRadiusCmd -> AvailableForStates(G4State_Idle);
	
  // Thickness of the second collimator
  ScollimatorThicknessCmd = new G4UIcmdWithADoubleAndUnit("/LaserDriven/EnergySelector/SecondCollimator/Thickness", this);
  ScollimatorThicknessCmd -> SetGuidance("Set the thickness of the second collimator");
  ScollimatorThicknessCmd -> SetParameterName("Size",false);
  ScollimatorThicknessCmd -> SetDefaultUnit("mm"); 
  ScollimatorThicknessCmd -> SetUnitCandidates("mm cm m");  
  ScollimatorThicknessCmd -> AvailableForStates(G4State_Idle);

  // Z Position of the second collimator hole 
  ScollimatorZpositionCmd = new G4UIcmdWithADoubleAndUnit("/LaserDriven/EnergySelector/SecondCollimator/zPosizion", this);
  ScollimatorZpositionCmd -> SetGuidance("Set the Z position of the second collimator hole as respect the internal vacuum chamber center axis");
  ScollimatorZpositionCmd -> SetParameterName("Size",false);
  ScollimatorZpositionCmd -> SetDefaultUnit("mm"); 
  ScollimatorZpositionCmd -> SetUnitCandidates("mm cm m");  
  ScollimatorZpositionCmd -> AvailableForStates(G4State_Idle);
		
  // THE SLIT
  //	
  // Thickness Slit
  SlitThicknessCmd = new G4UIcmdWithADoubleAndUnit("/LaserDriven/EnergySelector/Slit/thickness",this);
  SlitThicknessCmd -> SetGuidance("Set the X dimension of the Slit, the maximum value is 10 mm");
  SlitThicknessCmd -> SetParameterName("Size",false);
  SlitThicknessCmd -> SetDefaultUnit("mm"); 
  SlitThicknessCmd -> SetUnitCandidates("mm cm m");  
  SlitThicknessCmd -> AvailableForStates(G4State_Idle);
	
  //Hole dimension of the Slit (in Y direction)
  holeSlitDimensionYCmd = new G4UIcmdWithADoubleAndUnit("/LaserDriven/EnergySelector/Slit/HoleDimensionY",this);
  holeSlitDimensionYCmd -> SetGuidance("Set the Y dimension of the Slit Hole");
  holeSlitDimensionYCmd -> SetParameterName("Size",false);
  holeSlitDimensionYCmd -> SetDefaultUnit("mm"); 
  holeSlitDimensionYCmd -> SetUnitCandidates("mm cm m");  
  holeSlitDimensionYCmd -> AvailableForStates(G4State_Idle);
	
  // Hole dimension of the Slit (in Z direction)
  holeSlitDimensionZCmd = new G4UIcmdWithADoubleAndUnit("/LaserDriven/EnergySelector/Slit/HoleDimensionZ",this);
  holeSlitDimensionZCmd -> SetGuidance("Set the Z dimension of the external part of magnet 4");
  holeSlitDimensionZCmd -> SetParameterName("Size",false);
  holeSlitDimensionZCmd -> SetDefaultUnit("mm"); 
  holeSlitDimensionZCmd -> SetUnitCandidates("mm cm m");  
  holeSlitDimensionZCmd -> AvailableForStates(G4State_Idle);
	
  // Hole position of the Slit (in Z direction as respect the Slit body)
  slitHolePositionZCmd = new G4UIcmdWithADoubleAndUnit("/LaserDriven/EnergySelector/Slit/HolePositionZ", this);
  slitHolePositionZCmd -> SetGuidance("Set the Slit hole position in the Z direction as respect the Slit body center");
  slitHolePositionZCmd -> SetParameterName("Size",false);
  slitHolePositionZCmd -> SetDefaultUnit("mm"); 
  slitHolePositionZCmd -> SetUnitCandidates("mm cm m");  
  slitHolePositionZCmd -> AvailableForStates(G4State_Idle);

//  Quadrupole system DISABLE
   DisableQuadsCmd = new G4UIcmdWithoutParameter("/LaserDriven/Quadrupoles/DisableQuads", this);	
   DisableQuadsCmd -> SetGuidance("Disable the Quadrupole system.");
   DisableQuadsCmd -> SetGuidance("This command MUST be applied before \"beamOn\" ");
   DisableQuadsCmd -> AvailableForStates(G4State_Idle);			
}

/////////////////////////////////////////////////////////////////////////////
LaserDrivenBeamLineMessenger::~LaserDrivenBeamLineMessenger()
{ 

  delete laserDrivenDir;
  delete energySelectorDir;
  delete FcollimatorDir;
  delete ScollimatorDir;
  delete slitDir;
  delete quadrupoleDir;
  delete relativePosDir;
  delete DisableESSCmd;
  delete FcollimatorRadiusCmd;
  delete FcollimatorThicknessCmd;	
  delete FcollimatorZpositionCmd;
  delete ScollimatorRadiusCmd;
  delete ScollimatorThicknessCmd;	
  delete ScollimatorZpositionCmd;

  delete SlitThicknessCmd;
  delete holeSlitDimensionYCmd;
  delete holeSlitDimensionZCmd;
  delete slitHolePositionZCmd;

  delete DisableQuadsCmd;
		
}
/////////////////////////////////////////////////////////////////////////////
void LaserDrivenBeamLineMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if (command == DisableESSCmd)
    {
      laserDrivenMessengerPointer -> RemoveESS();
    }
  if( command == FcollimatorRadiusCmd )
    { 
      laserDrivenMessengerPointer -> SetFirstCollimatorRadius
	(FcollimatorRadiusCmd -> GetNewDoubleValue(newValue));
    }
  else if( command == FcollimatorThicknessCmd )
    { 
      laserDrivenMessengerPointer -> SetFirstCollimatorThickness
	(FcollimatorThicknessCmd -> GetNewDoubleValue(newValue));
    }
  else if( command == FcollimatorZpositionCmd )
    { 
      laserDrivenMessengerPointer -> SetFirstCollimatorPositionZ
	(FcollimatorZpositionCmd -> GetNewDoubleValue(newValue));
    }
  else if( command == ScollimatorRadiusCmd )
    { 
      laserDrivenMessengerPointer -> SetSecondCollimatorRadius
	(ScollimatorRadiusCmd -> GetNewDoubleValue(newValue));
    }
  else if( command == ScollimatorThicknessCmd )
    { 
      laserDrivenMessengerPointer -> SetSecondCollimatorThickness
	(ScollimatorThicknessCmd -> GetNewDoubleValue(newValue));
    }
  else if( command == ScollimatorZpositionCmd )
    { 
      laserDrivenMessengerPointer -> SetSecondCollimatorPositionZ
	(ScollimatorZpositionCmd -> GetNewDoubleValue(newValue));
    }
  else if( command == SlitThicknessCmd )
    { 
      laserDrivenMessengerPointer -> SetThicknessSlit
	(SlitThicknessCmd -> GetNewDoubleValue(newValue));
    }		
  else if( command == holeSlitDimensionYCmd ) 
    { 
      laserDrivenMessengerPointer -> SetSlitHoleDimensionY
	(holeSlitDimensionYCmd -> GetNewDoubleValue(newValue));
    }
  else if( command == holeSlitDimensionZCmd ) 
    { 
      laserDrivenMessengerPointer -> SetSlitHoleDimensionZ
	(holeSlitDimensionZCmd -> GetNewDoubleValue(newValue));
    }
  else if( command == slitHolePositionZCmd ) 
    { 
      laserDrivenMessengerPointer -> SetSlitHolePositionZ
	(slitHolePositionZCmd -> GetNewDoubleValue(newValue));
    }	
  else if (command==DisableQuadsCmd)
   {
      laserDrivenMessengerPointer -> RemoveQuads();
    }
  	
}
 
