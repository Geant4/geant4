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
/// \file medical/fanoCavity2/src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:fDetector(Det)
{ 
  fTestemDir = new G4UIdirectory("/testem/");
  fTestemDir->SetGuidance(" detector control.");
  
  fDetDir = new G4UIdirectory("/testem/det/");
  fDetDir->SetGuidance("detector construction commands");
      
  fWallMater = new G4UIcmdWithAString("/testem/det/wallMater",this);
  fWallMater->SetGuidance("Set material of the wall.");
  fWallMater->SetParameterName("wallMat",false);
  
  fWallThick = new G4UIcmdWithADoubleAndUnit("/testem/det/wallThickness",this);
  fWallThick->SetGuidance("Set tickness of the wall");
  fWallThick->SetParameterName("wallTick",false);
  fWallThick->SetRange("wallTick>0.");
  fWallThick->SetUnitCategory("Length");
  
  fCavThick = new G4UIcmdWithADoubleAndUnit("/testem/det/cavityThickness",this);
  fCavThick->SetGuidance("Set tickness of the cavity");
  fCavThick->SetParameterName("cavityTick",false);
  fCavThick->SetRange("cavityTick>0.");
  fCavThick->SetUnitCategory("Length");
  
  fWorldRadius = new G4UIcmdWithADoubleAndUnit("/testem/det/worldRadius",this);
  fWorldRadius->SetGuidance("Set radius of the cavity");
  fWorldRadius->SetParameterName("cavityRadius",false);
  fWorldRadius->SetRange("cavityRadius>0.");
  fWorldRadius->SetUnitCategory("Length");
        
  fUpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update",this);
  fUpdateCmd->SetGuidance("Update geometry.");
  fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fUpdateCmd->SetGuidance("if you changed geometrical value(s).");
  fUpdateCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fWallMater;
  
  delete fWallThick;  
  delete fCavThick;  
  delete fWorldRadius;
  
  delete fUpdateCmd;
  delete fDetDir;  
  delete fTestemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fWallMater )
   { fDetector->SetWallMaterial(newValue);}
   
  if( command == fWallThick )
   { fDetector->SetWallThickness(fWallThick->GetNewDoubleValue(newValue));}
   
  if( command == fCavThick )
   { fDetector->SetCavityThickness(fCavThick->GetNewDoubleValue(newValue));}
   
  if( command == fWorldRadius )
   { fDetector->SetWorldRadius(fWorldRadius->GetNewDoubleValue(newValue));}
                
  if( command == fUpdateCmd )
   { fDetector->UpdateGeometry();}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
