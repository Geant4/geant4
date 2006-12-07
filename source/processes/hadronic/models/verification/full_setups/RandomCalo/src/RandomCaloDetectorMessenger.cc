#include "RandomCaloDetectorMessenger.hh"

#include "RandomCaloDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "globals.hh"


RandomCaloDetectorMessenger::
RandomCaloDetectorMessenger( RandomCaloDetectorConstruction* myDet )
  : theDetector( myDet ) { 

  theDetectorDir = new G4UIdirectory( "/mydet/" );
  theDetectorDir->SetGuidance( "Detector control." );
  
  theFieldCommand = new G4UIcmdWithADoubleAndUnit( "/mydet/setField", this );  
  theFieldCommand->SetGuidance( "Define uniform magnetic field along Y." );
  theFieldCommand->SetGuidance( " -> in unit of  [Tesla]" );
  theFieldCommand->SetParameterName( "By", false );
  theFieldCommand->SetDefaultValue( 0.0 );
  theFieldCommand->SetUnitCategory( "Magnetic flux density" );
  theFieldCommand->AvailableForStates( G4State_PreInit, G4State_Idle );  
}


RandomCaloDetectorMessenger::~RandomCaloDetectorMessenger() {
  delete theFieldCommand;
  delete theDetectorDir;
}


void RandomCaloDetectorMessenger::
SetNewValue(G4UIcommand* command, G4String newValue) { 
  if ( command == theFieldCommand ) { 
    theDetector->SetMagField( theFieldCommand->GetNewDoubleValue(newValue) );
  }
}

