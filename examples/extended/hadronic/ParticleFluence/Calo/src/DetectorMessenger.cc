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
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger( DetectorConstruction* myDet ) : fDetector( myDet ) { 
  fDetectorDir = new G4UIdirectory( "/mydet/" );
  fDetectorDir->SetGuidance( "Detector control." );
  
  fFieldCommand = new G4UIcmdWithADoubleAndUnit( "/mydet/setField", this );  
  fFieldCommand->SetGuidance( "Define uniform magnetic field along Y." );
  fFieldCommand->SetGuidance( " -> in unit of  [Tesla]" );
  fFieldCommand->SetParameterName( "By", false );
  fFieldCommand->SetDefaultValue( 0.0 );
  fFieldCommand->SetUnitCategory( "Magnetic flux density" );
  fFieldCommand->AvailableForStates( G4State_PreInit, G4State_Idle );  

  fAbsorberMaterial = new G4UIcmdWithAString( "/mydet/absorberMaterial", this );
  fAbsorberMaterial->SetGuidance( "Choice of the absorber material:" );
  fAbsorberMaterial->SetGuidance( "   iron / copper / tungsten / lead / PbWO4 / uranium " );
  fAbsorberMaterial->SetParameterName( "choiceAbsorberMaterial", true );
  fAbsorberMaterial->SetDefaultValue( "iron" );
  fAbsorberMaterial->AvailableForStates( G4State_PreInit, G4State_Idle );

  fActiveMaterial = new G4UIcmdWithAString( "/mydet/activeMaterial", this );
  fActiveMaterial->SetGuidance( "Choice of the active material:" );
  fActiveMaterial->SetGuidance( "   scintillator / liquidArgon / PbWO4 / silicon / quartz " );
  fActiveMaterial->SetParameterName( "choiceActiveMaterial", true );
  fActiveMaterial->SetDefaultValue( "scintillator" );
  fActiveMaterial->AvailableForStates( G4State_PreInit, G4State_Idle );

  fIsCalHomogeneous = new G4UIcmdWithABool( "/mydet/isCalHomogeneous", this );
  fIsCalHomogeneous->SetParameterName( "choiceIsCalHomogeneous", true );
  fIsCalHomogeneous->SetGuidance( "Is the calorimeter homogeneous?" );
  fIsCalHomogeneous->SetGuidance( " -> yes|y|true|t|1 : Homogeneous calorimeter" );
  fIsCalHomogeneous->SetGuidance( " -> no|n|false|f|0 : Sampling calorimeter" );
  fIsCalHomogeneous->SetDefaultValue( false ); // default: sampling calorimeter
  fIsCalHomogeneous->AvailableForStates( G4State_PreInit, G4State_Idle );  
  
  fIsUnitInLambda = new G4UIcmdWithABool( "/mydet/isUnitInLambda", this );
  fIsUnitInLambda->SetParameterName( "choiceIsUnitInLambda", true );
  fIsUnitInLambda->SetGuidance( "Is unit for absorber length in lambda?" );
  fIsUnitInLambda->SetGuidance( " -> yes|y|true|t|1 : unit in lambda" );
  fIsUnitInLambda->SetGuidance( " -> no|n|false|f|0 : unit in [mm]" );
  fIsUnitInLambda->SetDefaultValue( false ); // default: unit in [mm].
  fIsUnitInLambda->AvailableForStates( G4State_PreInit, G4State_Idle );  
  
  fAbsorberTotalLength = new G4UIcmdWithADouble( "/mydet/absorberTotalLength", this );
  fAbsorberTotalLength->SetParameterName( "choiceAbsorberTotalLength", true );
  fAbsorberTotalLength->SetGuidance( "Absorber total length" );
  fAbsorberTotalLength->SetGuidance( " -> in unit of lambda or [mm]" );
  fAbsorberTotalLength->SetGuidance( " -> depending on value of choiceIsUnitInLambda" );
  fAbsorberTotalLength->SetDefaultValue( 2000.0 ); // default: 2 meters.
  fAbsorberTotalLength->AvailableForStates( G4State_PreInit, G4State_Idle );

  fCalorimeterRadius = new G4UIcmdWithADouble( "/mydet/calorimeterRadius", this );
  fCalorimeterRadius->SetParameterName( "choiceCalorimeterRadius", true );
  fCalorimeterRadius->SetGuidance( "Calorimeter Radius" );
  fCalorimeterRadius->SetGuidance( " -> in unit of lambda or [mm]" );
  fCalorimeterRadius->SetGuidance( " -> depending on value of choiceIsUnitInLambda" );
  fCalorimeterRadius->SetDefaultValue( 1000.0 ); // default: 1 meter.
  fCalorimeterRadius->AvailableForStates( G4State_PreInit, G4State_Idle );

  fActiveLayerNumber = new G4UIcmdWithAnInteger( "/mydet/activeLayerNumber", this );
  fActiveLayerNumber->SetParameterName( "choiceActiveLayerNumber", true );
  fActiveLayerNumber->SetGuidance( "Number of active layers" );
  fActiveLayerNumber->SetDefaultValue( 50 );
  fActiveLayerNumber->AvailableForStates( G4State_PreInit, G4State_Idle );

  fActiveLayerSize = new G4UIcmdWithADouble( "/mydet/activeLayerSize", this );
  fActiveLayerSize->SetParameterName( "choiceActiveLayerSize", true );
  fActiveLayerSize->SetGuidance( "Size (thickness) of the active layer, in [mm]" );
  fActiveLayerSize->SetDefaultValue( 4.0 ); // default: 4 millimeters.
  fActiveLayerSize->AvailableForStates( G4State_PreInit, G4State_Idle );

  fIsRadiusUnitInLambda = new G4UIcmdWithABool( "/mydet/isRadiusUnitInLambda", this );
  fIsRadiusUnitInLambda->SetParameterName( "choiceIsRadiusUnitInLambda", true );
  fIsRadiusUnitInLambda->SetGuidance( "Is unit of radius in lambda?" );
  fIsRadiusUnitInLambda->SetGuidance( " -> yes|y|true|t|1 : unit in lambda" );
  fIsRadiusUnitInLambda->SetGuidance( " -> no|n|false|f|0 : unit in [mm]" );
  fIsRadiusUnitInLambda->SetDefaultValue( false ); // default: unit in [mm].
  fIsRadiusUnitInLambda->AvailableForStates( G4State_PreInit, G4State_Idle );  
  
  fUpdateCommand = new G4UIcmdWithoutParameter( "/mydet/update", this);
  fUpdateCommand->SetGuidance( "Update calorimeter geometry." );
  fUpdateCommand->SetGuidance( "This command MUST be applied before \"beamOn\" " );
  fUpdateCommand->SetGuidance( "if you changed geometrical value(s)." );
  fUpdateCommand->AvailableForStates( G4State_Idle );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger() {
  delete fFieldCommand;
  delete fDetectorDir;
  delete fAbsorberMaterial;
  delete fActiveMaterial;
  delete fIsCalHomogeneous;
  delete fIsUnitInLambda;
  delete fAbsorberTotalLength;
  delete fCalorimeterRadius;
  delete fActiveLayerNumber;
  delete fActiveLayerSize;
  delete fIsRadiusUnitInLambda;
  delete fUpdateCommand;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue( G4UIcommand* command, G4String newValue ) { 
  if ( command == fFieldCommand ) { 
    fDetector->SetMagField( fFieldCommand->GetNewDoubleValue( newValue ) );
  }
  if ( command == fAbsorberMaterial ) { 
    fDetector->SetAbsorberMaterial( newValue );
  }
  if ( command == fActiveMaterial ) { 
    fDetector->SetActiveMaterial( newValue );
  }
  if ( command == fIsCalHomogeneous ) { 
    fDetector->SetIsCalHomogeneous( fIsCalHomogeneous->GetNewBoolValue( newValue ) );
  }
  if ( command == fIsUnitInLambda ) { 
    fDetector->SetIsUnitInLambda( fIsUnitInLambda->GetNewBoolValue( newValue ) );
  }
  if ( command == fAbsorberTotalLength ) { 
    fDetector->SetAbsorberTotalLength( fAbsorberTotalLength->GetNewDoubleValue( newValue ) );
  }
  if ( command == fCalorimeterRadius ) { 
    fDetector->SetCalorimeterRadius( fCalorimeterRadius->GetNewDoubleValue(newValue) );
  }
  if ( command == fActiveLayerNumber ) { 
    fDetector->SetActiveLayerNumber( fActiveLayerNumber->GetNewIntValue( newValue ) );
  }
  if ( command == fActiveLayerSize ) { 
    fDetector->SetActiveLayerSize( fActiveLayerSize->GetNewDoubleValue( newValue ) );
  }
  if ( command == fIsRadiusUnitInLambda ) { 
    fDetector->SetIsRadiusUnitInLambda( fIsRadiusUnitInLambda->GetNewBoolValue( newValue ) );
  }
  if ( command == fUpdateCommand ) {
    fDetector->UpdateGeometry();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
