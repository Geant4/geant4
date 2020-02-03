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
#include "G4RunManager.hh"

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
  sampleCmd->SetCandidates("Dolorite Anorthosite Mars1 IceBasalt");
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
  orbitHeightCmd->SetGuidance( "Set height of the spacecraft above Mercury Surface" );
  orbitHeightCmd->SetGuidance( "After this, /apparate/update must be executed before BeamOn" );
  orbitHeightCmd->SetGuidance( "Default: 400 km " );
  orbitHeightCmd->SetParameterName( "Spacecraft Altitude", true, true );
  orbitHeightCmd->SetDefaultUnit( "km" );
  orbitHeightCmd->SetUnitCategory( "Length" );
  orbitHeightCmd->AvailableForStates(G4State_Idle);
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
   { 
     //This triggers a full re-build of the geometry. The method in the 
     //geometry will take care of that.
     Detector->UpdateGeometry(); 
     return;
   }
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
 //Notify the run manager that the geometry has been modified
 G4RunManager::GetRunManager()->GeometryHasBeenModified();
 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....









