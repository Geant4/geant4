#include "globals.hh"
#include "G4VisAttributes.hh"
#include "MyDetectorConstruction.hh"
#include "G4Processor/GDMLProcessor.h"
#include "G4BooleanSolid.hh"
#include "G4CSGSolid.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "MyDetectorMessenger.hh"
#include "G4RunManager.hh"


MyDetectorConstruction::MyDetectorConstruction() :
  fieldMgr( 0 ) , uniformMagField( 0 ) , detectorMessenger( 0 ) {

  sxp.Initialize();
  config.SetURI( "cms.gdml" );                     //***LOOKHERE***
  config.SetSetupName( "Default" );
  sxp.Configure( &config );

  detectorMessenger = new MyDetectorMessenger( this );
}


MyDetectorConstruction::~MyDetectorConstruction() {
  delete uniformMagField;
  delete detectorMessenger;
  sxp.Finalize();
}


G4VPhysicalVolume* MyDetectorConstruction::Construct() { 

  fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  sxp.Run();
  
  fWorld = (G4VPhysicalVolume *) GDMLProcessor::GetInstance()->GetWorldVolume();

  fWorld->GetLogicalVolume()->SetVisAttributes (G4VisAttributes::Invisible);
  
  if( fWorld == 0 ) {
    G4Exception(
      "World volume not set properly check your setup selection criteria or GDML input!"
      );
  }

  return fWorld;
}


void MyDetectorConstruction::SetMagField( const G4double fieldValue ) {
  if ( uniformMagField ) {
    delete uniformMagField;
  }
  if ( std::abs( fieldValue ) > 0.0 ) {
    // Apply a global uniform magnetic field along the Z axis.
    // Notice that only if the magnetic field is not zero, the Geant4
    // transportion in field gets activated.

    uniformMagField = new G4UniformMagField( G4ThreeVector( 0.0, 0.0, fieldValue ) );

    fieldMgr->SetDetectorField( uniformMagField );
    fieldMgr->CreateChordFinder( uniformMagField );

    G4cout << G4endl
           << " *** SETTING MAGNETIC FIELD : fieldValue = " << fieldValue / tesla
           << " Tesla *** " << G4endl 
	   << G4endl;

  } 
}
