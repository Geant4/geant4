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

DetectorMessenger::DetectorMessenger( DetectorConstruction* myDet ) : G4UImessenger(),
                                                                      fDetector( myDet ) {
  fDetectorDir = new G4UIdirectory( "/mydet/" );
  fDetectorDir->SetGuidance( "Detector control." );

  fFieldCommand = new G4UIcmdWithADoubleAndUnit( "/mydet/setField", this );  
  fFieldCommand->SetGuidance( "Define uniform magnetic field along Y." );
  fFieldCommand->SetGuidance( " -> in unit of  [Tesla]" );
  fFieldCommand->SetParameterName( "By", false );
  fFieldCommand->SetDefaultValue( 0.0 );
  fFieldCommand->SetUnitCategory( "Magnetic flux density" );
  fFieldCommand->AvailableForStates( G4State_PreInit, G4State_Idle );  
  
  fTargetMaterial = new G4UIcmdWithAString( "/mydet/targetMaterial", this );
  fTargetMaterial->SetGuidance( "Choice of the target material:" );
  fTargetMaterial->SetGuidance( "   a Geant4 NIST material, e.g. G4_Be " );
  fTargetMaterial->SetParameterName( "choiceTargetMaterial", true );
  fTargetMaterial->SetDefaultValue( "G4_Be" );
  fTargetMaterial->AvailableForStates( G4State_PreInit, G4State_Idle );
  
  fTargetInnerRadius = new G4UIcmdWithADoubleAndUnit( "/mydet/targetInnerRadius", this );
  fTargetInnerRadius->SetParameterName( "choiceTargetInnerRadius", true );
  fTargetInnerRadius->SetGuidance( "Target inner radius" );
  fTargetInnerRadius->SetDefaultValue( 9.0 );
  fTargetInnerRadius->AvailableForStates( G4State_PreInit, G4State_Idle );

  fTargetOuterRadius = new G4UIcmdWithADoubleAndUnit( "/mydet/targetOuterRadius", this );
  fTargetOuterRadius->SetParameterName( "choiceTargetOuterRadius", true );
  fTargetOuterRadius->SetGuidance( "Target outer radius" );
  fTargetOuterRadius->SetDefaultValue( 11.0 );
  fTargetOuterRadius->AvailableForStates( G4State_PreInit, G4State_Idle );

  fUpdateCommand = new G4UIcmdWithoutParameter( "/mydet/update", this);
  fUpdateCommand->SetGuidance( "Update calorimeter geometry." );
  fUpdateCommand->SetGuidance( "This command MUST be applied before \"beamOn\" " );
  fUpdateCommand->SetGuidance( "if you changed geometrical value(s)." );
  fUpdateCommand->AvailableForStates( G4State_PreInit, G4State_Idle ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger() {
  delete fDetectorDir;
  delete fFieldCommand;
  delete fTargetMaterial;
  delete fTargetInnerRadius;
  delete fTargetOuterRadius;
  delete fUpdateCommand;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue( G4UIcommand* command, G4String newValue ) {
  if ( command == fFieldCommand ) { 
    fDetector->SetMagField( fFieldCommand->GetNewDoubleValue( newValue ) );
  }
  if ( command == fTargetMaterial ) { 
    fDetector->SetTargetMaterial( newValue );
  }
  if ( command == fTargetInnerRadius ) { 
    fDetector->SetTargetInnerRadius( fTargetInnerRadius->GetNewDoubleValue( newValue ) );
  }
  if ( command == fTargetOuterRadius ) { 
    fDetector->SetTargetOuterRadius( fTargetOuterRadius->GetNewDoubleValue( newValue ) );
  }
  if ( command == fUpdateCommand ) {
    fDetector->UpdateGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
