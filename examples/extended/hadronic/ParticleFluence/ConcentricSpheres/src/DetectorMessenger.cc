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

  fMaterialTracker = new G4UIcmdWithAString( "/mydet/material_tracker", this );
  fMaterialTracker->SetGuidance( "Choice of the material for the Tracker:" );
  fMaterialTracker->SetGuidance( "   a Geant4 NIST material, e.g. G4_Si " );
  fMaterialTracker->SetParameterName( "choiceMaterial", true );
  fMaterialTracker->SetDefaultValue( "G4_Si" );
  fMaterialTracker->AvailableForStates( G4State_PreInit, G4State_Idle );

  fMaterialEmCalo = new G4UIcmdWithAString( "/mydet/material_emCalo", this );
  fMaterialEmCalo->SetGuidance( "Choice of the material for the EM calo:" );
  fMaterialEmCalo->SetGuidance( "   a Geant4 NIST material, e.g. G4_Pb " );
  fMaterialEmCalo->SetParameterName( "choiceMaterial", true );
  fMaterialEmCalo->SetDefaultValue( "G4_Pb" );
  fMaterialEmCalo->AvailableForStates( G4State_PreInit, G4State_Idle );
  
  fMaterialHadCalo = new G4UIcmdWithAString( "/mydet/material_hadCalo", this );
  fMaterialHadCalo->SetGuidance( "Choice of the material for the HAD calo:" );
  fMaterialHadCalo->SetGuidance( "   a Geant4 NIST material, e.g. G4_Fe " );
  fMaterialHadCalo->SetParameterName( "choiceMaterial", true );
  fMaterialHadCalo->SetDefaultValue( "G4_Fe" );
  fMaterialHadCalo->AvailableForStates( G4State_PreInit, G4State_Idle );
  
  fInnerRadiusTracker = new G4UIcmdWithADoubleAndUnit( "/mydet/inner_radius_tracker", this );
  fInnerRadiusTracker->SetParameterName( "choiceInnerRadiusTracker", true );
  fInnerRadiusTracker->SetGuidance( "Inner radius of the Tracker" );
  fInnerRadiusTracker->SetDefaultValue( 100.0 );  // default: 10 cm.
  fInnerRadiusTracker->AvailableForStates( G4State_PreInit, G4State_Idle );

  fOuterRadiusTracker = new G4UIcmdWithADoubleAndUnit( "/mydet/outer_radius_tracker", this );
  fOuterRadiusTracker->SetParameterName( "choiceOuterRadiusTracker", true );
  fOuterRadiusTracker->SetGuidance( "Outer radius of the Tracker" );
  fOuterRadiusTracker->SetDefaultValue( 200.0 );  // default: 20 cm.
  fOuterRadiusTracker->AvailableForStates( G4State_PreInit, G4State_Idle );

  fInnerRadiusEmCalo = new G4UIcmdWithADoubleAndUnit( "/mydet/inner_radius_emCalo", this );
  fInnerRadiusEmCalo->SetParameterName( "choiceInnerRadiusEmCalo", true );
  fInnerRadiusEmCalo->SetGuidance( "Inner radius of the EM Calo" );
  fInnerRadiusEmCalo->SetDefaultValue( 300.0 );  // default: 30 cm.
  fInnerRadiusEmCalo->AvailableForStates( G4State_PreInit, G4State_Idle );

  fOuterRadiusEmCalo = new G4UIcmdWithADoubleAndUnit( "/mydet/outer_radius_emCalo", this );
  fOuterRadiusEmCalo->SetParameterName( "choiceOuterRadiusEmCalo", true );
  fOuterRadiusEmCalo->SetGuidance( "Outer radius of the EM Calo" );
  fOuterRadiusEmCalo->SetDefaultValue( 600.0 );  // default: 60 cm.
  fOuterRadiusEmCalo->AvailableForStates( G4State_PreInit, G4State_Idle );
  
  fInnerRadiusHadCalo = new G4UIcmdWithADoubleAndUnit( "/mydet/inner_radius_hadCalo", this );
  fInnerRadiusHadCalo->SetParameterName( "choiceInnerRadiusHadCalo", true );
  fInnerRadiusHadCalo->SetGuidance( "Inner radius of the HAD Calo" );
  fInnerRadiusHadCalo->SetDefaultValue( 700.0 );  // default: 70 cm.
  fInnerRadiusHadCalo->AvailableForStates( G4State_PreInit, G4State_Idle );

  fOuterRadiusHadCalo = new G4UIcmdWithADoubleAndUnit( "/mydet/outer_radius_hadCalo", this );
  fOuterRadiusHadCalo->SetParameterName( "choiceOuterRadiusHadCalo", true );
  fOuterRadiusHadCalo->SetGuidance( "Outer radius of the HAD Calo" );
  fOuterRadiusHadCalo->SetDefaultValue( 1700.0 );  // default: 170 cm.
  fOuterRadiusHadCalo->AvailableForStates( G4State_PreInit, G4State_Idle );
  
  fUpdateCommand = new G4UIcmdWithoutParameter( "/mydet/update", this);
  fUpdateCommand->SetGuidance( "Update geometry." );
  fUpdateCommand->SetGuidance( "This command MUST be applied before \"beamOn\" " );
  fUpdateCommand->SetGuidance( "if you changed geometrical value(s)." );
  fUpdateCommand->AvailableForStates( G4State_Idle ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger() {
  delete fDetectorDir;
  delete fMaterialTracker;
  delete fMaterialEmCalo;
  delete fMaterialHadCalo;
  delete fInnerRadiusTracker;
  delete fOuterRadiusTracker;
  delete fInnerRadiusEmCalo;
  delete fOuterRadiusEmCalo;
  delete fInnerRadiusHadCalo;
  delete fOuterRadiusHadCalo;
  delete fUpdateCommand;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue( G4UIcommand* command, G4String newValue ) {
  if ( command == fMaterialTracker ) { 
    fDetector->SetMaterialTracker( newValue );
  }
  if ( command == fMaterialEmCalo ) { 
    fDetector->SetMaterialEmCalo( newValue );
  }
  if ( command == fMaterialHadCalo ) { 
    fDetector->SetMaterialHadCalo( newValue );
  }  
  if ( command == fInnerRadiusTracker ) {
    fDetector->SetInnerRadiusTracker( fInnerRadiusTracker->GetNewDoubleValue( newValue ) );
  }
  if ( command == fOuterRadiusTracker ) {
    fDetector->SetOuterRadiusTracker( fOuterRadiusTracker->GetNewDoubleValue( newValue ) );
  }
  if ( command == fInnerRadiusEmCalo ) {
    fDetector->SetInnerRadiusEmCalo( fInnerRadiusEmCalo->GetNewDoubleValue( newValue ) );
  }
  if ( command == fOuterRadiusEmCalo ) {
    fDetector->SetOuterRadiusEmCalo( fOuterRadiusEmCalo->GetNewDoubleValue( newValue ) );
  }
  if ( command == fInnerRadiusHadCalo ) {
    fDetector->SetInnerRadiusHadCalo( fInnerRadiusHadCalo->GetNewDoubleValue( newValue ) );
  }
  if ( command == fOuterRadiusHadCalo ) {
    fDetector->SetOuterRadiusHadCalo( fOuterRadiusHadCalo->GetNewDoubleValue( newValue ) );
  }
  if ( command == fUpdateCommand ) {
    fDetector->UpdateGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
