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
  fWallMater->AvailableForStates(G4State_PreInit);
  
  fWallThick = new G4UIcmdWithADoubleAndUnit("/testem/det/wallThickness",this);
  fWallThick->SetGuidance("Set tickness of the wall");
  fWallThick->SetParameterName("wallTick",false);
  fWallThick->SetRange("wallTick>0.");
  fWallThick->SetUnitCategory("Length");
  fWallThick->AvailableForStates(G4State_PreInit);
    
  fCavMater = new G4UIcmdWithAString("/testem/det/cavityMater",this);
  fCavMater->SetGuidance("Set material of the cavity.");
  fCavMater->SetParameterName("cavMat",false);
  fCavMater->AvailableForStates(G4State_PreInit);
  
  fCavThick = new G4UIcmdWithADoubleAndUnit("/testem/det/cavityThickness",this);
  fCavThick->SetGuidance("Set tickness of the cavity");
  fCavThick->SetParameterName("cavityTick",false);
  fCavThick->SetRange("cavityTick>0.");
  fCavThick->SetUnitCategory("Length");
  fCavThick->AvailableForStates(G4State_PreInit);
  
  fCavRadius = new G4UIcmdWithADoubleAndUnit("/testem/det/cavityRadius",this);
  fCavRadius->SetGuidance("Set radius of the cavity");
  fCavRadius->SetParameterName("cavityRadius",false);
  fCavRadius->SetRange("cavityRadius>0.");
  fCavRadius->SetUnitCategory("Length");
  fCavRadius->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fWallMater;
  delete fWallThick;
  
  delete fCavMater;
  delete fCavThick;
  delete fCavRadius;
  
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
      
  if( command == fCavMater )
   { fDetector->SetCavityMaterial(newValue);}
   
  if( command == fCavThick )
   { fDetector->SetCavityThickness(fCavThick->GetNewDoubleValue(newValue));}
   
  if( command == fCavRadius )
   { fDetector->SetCavityRadius(fCavRadius->GetNewDoubleValue(newValue));}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
