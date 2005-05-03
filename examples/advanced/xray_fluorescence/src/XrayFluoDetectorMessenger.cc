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
// $Id: XrayFluoDetectorMessenger.cc
// GEANT4 tag $Name: xray_fluo-V03-02-00
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#include "XrayFluoDetectorMessenger.hh"
#include "XrayFluoDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoDetectorMessenger::XrayFluoDetectorMessenger(XrayFluoDetectorConstruction * Det)
:Detector(Det)
{ 
  detDir = new G4UIdirectory("/apparate/");
  detDir->SetGuidance("detector control.");

  UpdateCmd = new G4UIcmdWithoutParameter("/apparate/update",this);
  UpdateCmd->SetGuidance("Update apparate geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s): /apparate/GrainDiameter and /apparate/sampleGranularity");

  UpdateCmd->AvailableForStates(G4State_Idle);

  sampleCmd = new G4UIcmdWithAString("/apparate/sampleMaterial",this);
  sampleCmd->SetGuidance("select a diferent material for the sample");
  sampleCmd->SetParameterName("material",true);
  sampleCmd->SetDefaultValue("mars1");
  sampleCmd->SetCandidates("Dolorite Iron Silicon Aluinium Oxigen Titanium Tin Lead Neodimium Magnesium Copper Anorthosite Mars1 IceBasalt Silver Gold Caesium Potassium Manganese Phosphorus Sulphur Calcium Sodium Uranium HPGe");
  sampleCmd->AvailableForStates(G4State_Idle);

  detectorCmd = new G4UIcmdWithAString("/apparate/detector",this);
  detectorCmd->SetGuidance("select a diferent detectorType");
  detectorCmd->SetParameterName("detector",true);
  detectorCmd->SetDefaultValue("sili");
  detectorCmd->SetCandidates("sili hpge");
  detectorCmd->AvailableForStates(G4State_Idle);
  
  grainDiaCmd = new G4UIcmdWithADoubleAndUnit( "/apparate/GrainDiameter",this );
  grainDiaCmd->SetGuidance( "Set diameter of grains" );
  grainDiaCmd->SetGuidance( "After this, /apparate/update must be executed before BeamOn" );
  grainDiaCmd->SetGuidance( "Default: 0.5 mm " );
  grainDiaCmd->SetParameterName( "Grain Diameter", true, true );
  grainDiaCmd->SetDefaultUnit( "mm" );
  grainDiaCmd->SetUnitCategory( "Length" );
  grainDiaCmd->AvailableForStates(G4State_Idle);

  granularityFlagCmd= new G4UIcmdWithABool("/apparate/sampleGranularity",this);
  granularityFlagCmd->SetGuidance("Set if sample granularity is present");
  granularityFlagCmd->SetGuidance( "After this, /apparate/update must be executed before BeamOn" );
  granularityFlagCmd->SetParameterName("Granularity Flag",true);
  granularityFlagCmd->SetDefaultValue(false);
  granularityFlagCmd->AvailableForStates(G4State_Idle);

  OhmicPosThicknessCmd = new G4UIcmdWithADoubleAndUnit( "/apparate/ohmicPosThickness",this );
  OhmicPosThicknessCmd->SetGuidance( "Changes thickness of the detector anode" );
  OhmicPosThicknessCmd->SetGuidance( "After this, /apparate/update must be executed before BeamOn" );
  OhmicPosThicknessCmd->SetGuidance( "Default: 0.005 mm " );
  OhmicPosThicknessCmd->SetParameterName( "OhmicPos Thickness", true, true );
  OhmicPosThicknessCmd->SetDefaultUnit( "mm" );
  OhmicPosThicknessCmd->SetUnitCategory( "Length" );
  OhmicPosThicknessCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

 XrayFluoDetectorMessenger::~XrayFluoDetectorMessenger()
{
  delete UpdateCmd;
  delete detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
 if( command == UpdateCmd )
   { Detector->UpdateGeometry(); }

 else if ( command == sampleCmd )
   { Detector->SetSampleMaterial(newValue);}

 else if ( command == detectorCmd )
   { Detector->SetDetectorType(newValue);}

 else if ( command == grainDiaCmd )
   {  
     G4double newSize = grainDiaCmd->GetNewDoubleValue(newValue);  
     Detector->SetGrainDia(newSize);
   }

 else if ( command == granularityFlagCmd )
   { 
     Detector->DeleteGrainObjects();
     G4bool newGranFlag = granularityFlagCmd->GetNewBoolValue(newValue);
     Detector->SetSampleGranularity(newGranFlag);
   }
 else if ( command == OhmicPosThicknessCmd )
   {  
     G4double newSize = OhmicPosThicknessCmd->GetNewDoubleValue(newValue);  
     Detector->SetOhmicPosThickness(newSize);
   }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....









