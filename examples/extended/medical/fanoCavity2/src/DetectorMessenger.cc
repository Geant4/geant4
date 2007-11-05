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
// $Id: DetectorMessenger.cc,v 1.2 2007-11-05 13:19:16 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
:Detector(Det)
{ 
  testemDir = new G4UIdirectory("/testem/");
  testemDir->SetGuidance(" detector control.");
  
  detDir = new G4UIdirectory("/testem/det/");
  detDir->SetGuidance("detector construction commands");
      
  wallMater = new G4UIcmdWithAString("/testem/det/wallMater",this);
  wallMater->SetGuidance("Set material of the wall.");
  wallMater->SetParameterName("wallMat",false);
  
  wallThick = new G4UIcmdWithADoubleAndUnit("/testem/det/wallThickness",this);
  wallThick->SetGuidance("Set tickness of the wall");
  wallThick->SetParameterName("wallTick",false);
  wallThick->SetRange("wallTick>0.");
  wallThick->SetUnitCategory("Length");
  
  cavThick = new G4UIcmdWithADoubleAndUnit("/testem/det/cavityThickness",this);
  cavThick->SetGuidance("Set tickness of the cavity");
  cavThick->SetParameterName("cavityTick",false);
  cavThick->SetRange("cavityTick>0.");
  cavThick->SetUnitCategory("Length");
  
  worldRadius = new G4UIcmdWithADoubleAndUnit("/testem/det/worldRadius",this);
  worldRadius->SetGuidance("Set radius of the cavity");
  worldRadius->SetParameterName("cavityRadius",false);
  worldRadius->SetRange("cavityRadius>0.");
  worldRadius->SetUnitCategory("Length");
        
  UpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update",this);
  UpdateCmd->SetGuidance("Update geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete wallMater;
  
  delete wallThick;  
  delete cavThick;  
  delete worldRadius;
  
  delete UpdateCmd;
  delete detDir;  
  delete testemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == wallMater )
   { Detector->SetWallMaterial(newValue);}
   
  if( command == wallThick )
   { Detector->SetWallThickness(wallThick->GetNewDoubleValue(newValue));}
   
  if( command == cavThick )
   { Detector->SetCavityThickness(cavThick->GetNewDoubleValue(newValue));}
   
  if( command == worldRadius )
   { Detector->SetWorldRadius(worldRadius->GetNewDoubleValue(newValue));}
                
  if( command == UpdateCmd )
   { Detector->UpdateGeometry();}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
