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
// $Id: XrayFluoMercuryDetectorMessenger.cc
// GEANT4 tag $Name: 
//
// Author: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
//
//  18 Sep 2003  Alfonso Mantero created
//
// -------------------------------------------------------------------


#include "XrayFluoMercuryDetectorMessenger.hh"
#include "XrayFluoMercuryDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoMercuryDetectorMessenger::XrayFluoMercuryDetectorMessenger(XrayFluoMercuryDetectorConstruction * Det)
:Detector(Det)
{ 
  detDir = new G4UIdirectory("/apparate/");
  detDir->SetGuidance("detector control.");

  UpdateCmd = new G4UIcmdWithoutParameter("/apparate/update",this);
  UpdateCmd->SetGuidance("Update apparate geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s): /apparate/GrainDiameter and /apparate/sampleGranularity");

  UpdateCmd->AvailableForStates(G4State_Idle);

  sampleCmd = new G4UIcmdWithAString("/apparate/mercuryMaterial",this);
  sampleCmd->SetGuidance("select a diferent material for the mercury");
  sampleCmd->SetParameterName("material",true);
  sampleCmd->SetDefaultValue("mars1");
  sampleCmd->SetCandidates("Dolorite Anorthosite Mars1");
  sampleCmd->AvailableForStates(G4State_Idle);

  detectorCmd = new G4UIcmdWithAString("/apparate/detector",this);
  detectorCmd->SetGuidance("select a diferent detectorType");
  detectorCmd->SetParameterName("detector",true);
  detectorCmd->SetDefaultValue("sili");
  detectorCmd->SetCandidates("sili hpge");
  detectorCmd->AvailableForStates(G4State_Idle);
  
  latitudeAngleCmd = new G4UIcmdWithADoubleAndUnit( "/apparate/latitude",this );
  latitudeAngleCmd->SetGuidance( "Set latitude angle of the spacecraft" );
  latitudeAngleCmd->SetGuidance( "After this, /apparate/update must be executed before BeamOn" );
  latitudeAngleCmd->SetGuidance( "Default: 45 deg " );
  latitudeAngleCmd->SetParameterName( "Latitude Angle", true, true );
  latitudeAngleCmd->SetDefaultUnit( "deg" );
  latitudeAngleCmd->SetUnitCategory( "Angle" );
  latitudeAngleCmd->AvailableForStates(G4State_Idle);

  orbitHeightCmd = new G4UIcmdWithADoubleAndUnit( "/apparate/orbitHeight",this );
  orbitHeightCmd->SetGuidance( "Set height of the spacecraft above Mercuey Surface" );
  orbitHeightCmd->SetGuidance( "After this, /apparate/update must be executed before BeamOn" );
  orbitHeightCmd->SetGuidance( "Default: 400 km " );
  orbitHeightCmd->SetParameterName( "Spacecraft Altitude", true, true );
  orbitHeightCmd->SetDefaultUnit( "km" );
  orbitHeightCmd->SetUnitCategory( "Length" );
  orbitHeightCmd->AvailableForStates(G4State_Idle);



//   granularityFlagCmd= new G4UIcmdWithABool("/apparate/sampleGranularity",this);
//   granularityFlagCmd->SetGuidance("Set if sample granularity is present");
//   granularityFlagCmd->SetGuidance( "After this, /apparate/update must be executed before BeamOn" );
//   granularityFlagCmd->SetParameterName("Granularity Flag",true);
//   granularityFlagCmd->SetDefaultValue(false);
//   granularityFlagCmd->AvailableForStates(G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

 XrayFluoMercuryDetectorMessenger::~XrayFluoMercuryDetectorMessenger()
{
  delete UpdateCmd;
  delete detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoMercuryDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
 if( command == UpdateCmd )
   { Detector->UpdateGeometry(); }

 else if ( command == sampleCmd )
   { Detector->SetMercuryMaterial(newValue);}

 else if ( command == detectorCmd )
   { Detector->SetDetectorType(newValue);}

 else if ( command == latitudeAngleCmd )
   {  
     G4double newAngle = latitudeAngleCmd->GetNewDoubleValue(newValue);  
     Detector->SetLatitude(newAngle);
   }

 else if ( command == orbitHeightCmd )
   {  
     G4double newAngle = orbitHeightCmd->GetNewDoubleValue(newValue);  
     Detector->SetOribitHeight(newAngle);
   }


//  else if ( command == granularityFlagCmd )
//    { 
//      Detector->DeleteGrainObjects();
//      G4bool newGranFlag = granularityFlagCmd->GetNewBoolValue(newValue);
//      Detector->SetMercuryGranularity(newGranFlag);
//    }
 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....









