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
//  29 Aug 2003  Alfonso Mantero created
//
// -------------------------------------------------------------------


#include "XrayFluoPlaneDetectorMessenger.hh"
#include "XrayFluoPlaneDetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPlaneDetectorMessenger::XrayFluoPlaneDetectorMessenger(XrayFluoPlaneDetectorConstruction * Det)
:Detector(Det)
{ 
  detDir = new G4UIdirectory("/apparate/");
  detDir->SetGuidance("detector control.");

  UpdateCmd = new G4UIcmdWithoutParameter("/apparate/update",this);
  UpdateCmd->SetGuidance("Update apparate geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s): /apparate/GrainDiameter and /apparate/sampleGranularity");

  UpdateCmd->AvailableForStates(G4State_Idle);

  sampleCmd = new G4UIcmdWithAString("/apparate/planeMaterial",this);
  sampleCmd->SetGuidance("select a diferent material for the plane");
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

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

 XrayFluoPlaneDetectorMessenger::~XrayFluoPlaneDetectorMessenger()
{
  delete UpdateCmd;
  delete detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPlaneDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
 if( command == UpdateCmd )
   { 
     //This triggers a full re-build of the geometry. The method in the 
     //geometry will take care of that.
     Detector->UpdateGeometry(); 
     return;
   }
 else if ( command == sampleCmd )
   { 
     Detector->SetPlaneMaterial(newValue);
   }
 else if ( command == detectorCmd )
   { 
     Detector->SetDetectorType(newValue);
   }
 else if ( command == grainDiaCmd )
   {  
     G4double newSize = grainDiaCmd->GetNewDoubleValue(newValue);  
     Detector->SetGrainDia(newSize);
   }
 else if ( command == granularityFlagCmd )
   { 
     Detector->DeleteGrainObjects();
     G4bool newGranFlag = granularityFlagCmd->GetNewBoolValue(newValue);
     Detector->SetPlaneGranularity(newGranFlag);
   }
 //Notify the run manager that the geometry has been modified
 G4RunManager::GetRunManager()->GeometryHasBeenModified();
 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....









