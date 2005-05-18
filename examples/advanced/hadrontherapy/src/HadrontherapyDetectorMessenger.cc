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
// $Id: HadrontherapyDetectorMessenger.cc,v 2.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G. Candiano, G.A.P. Cirrone, F. Di Rosa, G. Russo 
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------
#include "HadrontherapyDetectorMessenger.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"

HadrontherapyDetectorMessenger::HadrontherapyDetectorMessenger(
                                           HadrontherapyDetectorConstruction* HadrontherapyDet)
:HadrontherapyDetector(HadrontherapyDet)
{ 
  modulatorDir = new G4UIdirectory("/modulator/");
  modulatorDir -> SetGuidance("Command to rotate the modulator wheel");
  
  beamLineDir = new G4UIdirectory("/beamLine/");
  beamLineDir->SetGuidance("set specification of range shifter");  

  RangeShifterDir = new G4UIdirectory("/beamLine/RangeShifter/");
  RangeShifterDir->SetGuidance("set specification of range shifter");  

  FirstScatteringFoilDir = new G4UIdirectory("/beamLine/ScatteringFoil1/");
  FirstScatteringFoilDir->SetGuidance("set specification of first scattering foil");  
 
  SecondScatteringFoilDir = new G4UIdirectory("/beamLine/ScatteringFoil2/");
  SecondScatteringFoilDir->SetGuidance("set specification of second scattering foil");  
 
  StopperDir = new G4UIdirectory("/beamLine/Stopper/");
  StopperDir->SetGuidance("set specification of stopper");  

  FinalCollimatorDir = new G4UIdirectory("/beamLine/FinalCollimator/");
  FinalCollimatorDir->SetGuidance("set specification of final collimator");  

   ModulatorAngleCmd = new G4UIcmdWithADoubleAndUnit("/modulator/angle",this);
   ModulatorAngleCmd->SetGuidance("Set Modulator Angle");
   ModulatorAngleCmd->SetParameterName("Size",false);
   ModulatorAngleCmd->SetRange("Size>=0.");
   ModulatorAngleCmd->SetUnitCategory("Angle");  
   ModulatorAngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
   RangeShifterMatCmd = new G4UIcmdWithAString("/beamLine/RangeShifter/RSMat",this);
   RangeShifterMatCmd->SetGuidance("Set material of range shifter");
   RangeShifterMatCmd->SetParameterName("choice",false);
   RangeShifterMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
   RangeShifterBox_xCmd = new G4UIcmdWithADoubleAndUnit("/beamLine/RangeShifter/thickness",this);
   RangeShifterBox_xCmd->SetGuidance("Set half of the thickness of range shifter");
   RangeShifterBox_xCmd->SetParameterName("Size",false);
   RangeShifterBox_xCmd->SetDefaultUnit("mm");  
   RangeShifterBox_xCmd->SetUnitCandidates("mm cm m");  
   RangeShifterBox_xCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
   RangeShifterBoxPosition_xCmd = new G4UIcmdWithADoubleAndUnit("/beamLine/RangeShifter/position",this);
   RangeShifterBoxPosition_xCmd->SetGuidance("Set position of range shifter");
   RangeShifterBoxPosition_xCmd->SetParameterName("Size",false);
   RangeShifterBoxPosition_xCmd->SetDefaultUnit("mm");  
   RangeShifterBoxPosition_xCmd->SetUnitCandidates("mm cm m");  
   RangeShifterBoxPosition_xCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
   FirstScatteringFoil_xCmd = new G4UIcmdWithADoubleAndUnit("/beamLine/ScatteringFoil1/thickness",this);
   FirstScatteringFoil_xCmd->SetGuidance("Set thickness of first scattering foil");
   FirstScatteringFoil_xCmd->SetParameterName("Size",false);
   FirstScatteringFoil_xCmd->SetDefaultUnit("mm");  
   FirstScatteringFoil_xCmd->SetUnitCandidates("mm cm m");  
   FirstScatteringFoil_xCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
   SecondScatteringFoil_xCmd = new G4UIcmdWithADoubleAndUnit("/beamLine/ScatteringFoil2/thickness",this);
   SecondScatteringFoil_xCmd->SetGuidance("Set thickness of second scattering foil");
   SecondScatteringFoil_xCmd->SetParameterName("Size",false);
   SecondScatteringFoil_xCmd->SetDefaultUnit("mm");  
   SecondScatteringFoil_xCmd->SetUnitCandidates("mm cm m");  
   SecondScatteringFoil_xCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
   outerRadiusStopperCmd = new G4UIcmdWithADoubleAndUnit("/beamLine/Stopper/outRadius",this);
   outerRadiusStopperCmd->SetGuidance("Set size of outer radius");
   outerRadiusStopperCmd->SetParameterName("Size",false);
   outerRadiusStopperCmd->SetDefaultUnit("mm");  
   outerRadiusStopperCmd->SetUnitCandidates("mm cm m");  
   outerRadiusStopperCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
   innerRadiusFinalCollimatorCmd = new G4UIcmdWithADoubleAndUnit("/beamLine/FinalCollimator/halfInnerRad",this);
   innerRadiusFinalCollimatorCmd->SetGuidance("Set size of inner radius ( max 21.5 mm)");
   innerRadiusFinalCollimatorCmd->SetParameterName("Size",false);
   innerRadiusFinalCollimatorCmd->SetDefaultUnit("mm");  
   innerRadiusFinalCollimatorCmd->SetUnitCandidates("mm cm m");  
   innerRadiusFinalCollimatorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

HadrontherapyDetectorMessenger::~HadrontherapyDetectorMessenger()
{
  delete ModulatorAngleCmd;
  delete detDir;
  delete modulatorDir;  
  delete beamLineDir;  
  delete RangeShifterDir;  
  delete FirstScatteringFoilDir;  
  delete SecondScatteringFoilDir; 
  delete StopperDir; 
  delete FinalCollimatorDir; 
  delete RangeShifterMatCmd;  
  delete RangeShifterBox_xCmd;  
  delete RangeShifterBoxPosition_xCmd;  
  delete FirstScatteringFoil_xCmd;  
  delete SecondScatteringFoil_xCmd;  
  delete outerRadiusStopperCmd;  
  delete innerRadiusFinalCollimatorCmd;  



}

void HadrontherapyDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  if( command == ModulatorAngleCmd )
  { HadrontherapyDetector->SetModulatorAngle(ModulatorAngleCmd->GetNewDoubleValue(newValue));}


  if( command == RangeShifterMatCmd )
    { HadrontherapyDetector->SetRSMaterial(newValue);}

 if( command == RangeShifterBox_xCmd )
  { HadrontherapyDetector->SetRangeShifterX(RangeShifterBox_xCmd->GetNewDoubleValue(newValue));}

 if( command == RangeShifterBoxPosition_xCmd )
  { HadrontherapyDetector->SetRangeShifterXposition(RangeShifterBoxPosition_xCmd->GetNewDoubleValue(newValue));}


 if( command == FirstScatteringFoil_xCmd )
  { HadrontherapyDetector->SetFirstScatteringFoil(FirstScatteringFoil_xCmd->GetNewDoubleValue(newValue));}

 if( command == SecondScatteringFoil_xCmd )
  { HadrontherapyDetector->SetSecondScatteringFoil(SecondScatteringFoil_xCmd->GetNewDoubleValue(newValue));}

 if( command == outerRadiusStopperCmd )
  { HadrontherapyDetector->SetOuterRadiusStopper(outerRadiusStopperCmd->GetNewDoubleValue(newValue));}


 if( command == innerRadiusFinalCollimatorCmd )
  { HadrontherapyDetector->SetInnerRadiusFinalCollimator(innerRadiusFinalCollimatorCmd->GetNewDoubleValue(newValue));}

}

