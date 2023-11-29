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
/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger( DetectorConstruction* myDet ) : fDetector( myDet ) {
  fDetectorDir = new G4UIdirectory( "/mydet/" );
  fDetectorDir->SetGuidance( "Detector control." );

  fMaterial = new G4UIcmdWithAString( "/mydet/material", this );
  fMaterial->SetGuidance( "Choice of the material:" );
  fMaterial->SetGuidance( "   a Geant4 NIST material, e.g. G4_Fe " );
  fMaterial->SetParameterName( "choiceMaterial", true );
  fMaterial->SetDefaultValue( "G4_Fe" );
  fMaterial->AvailableForStates( G4State_PreInit, G4State_Idle );
  
  fThickness = new G4UIcmdWithADoubleAndUnit( "/mydet/thickness", this );
  fThickness->SetParameterName( "choiceThickness", true );
  fThickness->SetGuidance( "Target thickness" );
  fThickness->SetDefaultValue( 1000.0 ); // default: 1 meter.
  fThickness->AvailableForStates( G4State_PreInit, G4State_Idle );

  fDiameter = new G4UIcmdWithADoubleAndUnit( "/mydet/diameter", this );
  fDiameter->SetParameterName( "choiceDiameter", true );
  fDiameter->SetGuidance( "Target diameter" );
  fDiameter->SetDefaultValue( 1000.0 ); // default: 1 meter.
  fDiameter->AvailableForStates( G4State_PreInit, G4State_Idle );

  fUpdateCommand = new G4UIcmdWithoutParameter( "/mydet/update", this);
  fUpdateCommand->SetGuidance( "Update calorimeter geometry." );
  fUpdateCommand->SetGuidance( "This command MUST be applied before \"beamOn\" " );
  fUpdateCommand->SetGuidance( "if you changed geometrical value(s)." );
  fUpdateCommand->AvailableForStates( G4State_Idle ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger() {
  delete fDetectorDir;
  delete fMaterial;
  delete fThickness;
  delete fDiameter;
  delete fUpdateCommand;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue( G4UIcommand* command, G4String newValue ) {
  if ( command == fMaterial ) { 
    fDetector->SetMaterial( newValue );
  }
  if ( command == fThickness ) { 
    fDetector->SetThickness( fThickness->GetNewDoubleValue( newValue ) );
  }
  if ( command == fDiameter ) { 
    fDetector->SetDiameter( fDiameter->GetNewDoubleValue(newValue) );
  }
  if ( command == fUpdateCommand ) {
    fDetector->UpdateGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
