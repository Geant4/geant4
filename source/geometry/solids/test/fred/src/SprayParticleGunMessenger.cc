//
// SprayParticleGunMessenger.cc
//
// Implementation of SprayParticleGun's UI interface
//

#include "SprayParticleGunMessenger.hh"
#include "SprayParticleGun.hh"
#include "G4ThreeVector.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4ios.hh"

//
// Constructor
//
SprayParticleGunMessenger::SprayParticleGunMessenger( SprayParticleGun *ourGun )
{
	gun = ourGun;

	gunDirectory = new G4UIdirectory( "/sprayGun/" );
	gunDirectory->SetGuidance( "Spray Particle Gun control commands." );
	
	positionCmd = new G4UIcmdWith3VectorAndUnit( "/sprayGun/position", this );
	positionCmd->SetGuidance( "Set starting position of particle(s)." );
	positionCmd->SetParameterName( "X", "Y", "Z", true, true );
	positionCmd->SetDefaultUnit( "cm" );
	
	xSprayCmd = new G4UIcmdWithAnInteger( "/sprayGun/xSpray", this );
	xSprayCmd->SetGuidance( "X spray parameter" );
	xSprayCmd->SetParameterName( "xSpray", true, true );
//	xSprayCmd->SetRange( "N<8" );
	
	ySprayCmd = new G4UIcmdWithAnInteger( "/sprayGun/ySpray", this );
	ySprayCmd->SetGuidance( "Y spray parameter" );
	ySprayCmd->SetParameterName( "ySpray", true, true );
//	ySprayCmd->SetRange( "N<8" );
	
	zSprayCmd = new G4UIcmdWithAnInteger( "/sprayGun/zSpray", this );
	zSprayCmd->SetGuidance( "Z spray parameter" );
	zSprayCmd->SetParameterName( "zSpray", true, true );
//	zSprayCmd->SetRange( "N<8" );
}

//
// Destructor
//
SprayParticleGunMessenger::~SprayParticleGunMessenger()
{
	// is the order important?
	delete positionCmd;
	delete xSprayCmd;
	delete ySprayCmd;
	delete zSprayCmd;
	delete gunDirectory;
}


//
// SetNewValue
//
// Called by the UI when user requests a change
//
void SprayParticleGunMessenger::SetNewValue( G4UIcommand *command, G4String newValues )
{
	if (command == positionCmd) {
		gun->SetPosition( positionCmd->GetNew3VectorValue( newValues ) );
	}
	else if (command == xSprayCmd) {
		gun->SetXSpray( xSprayCmd->GetNewIntValue( newValues ) );
	}
	else if (command == ySprayCmd) {
		gun->SetYSpray( ySprayCmd->GetNewIntValue( newValues ) );
	}
	else if (command == zSprayCmd) {
		gun->SetZSpray( zSprayCmd->GetNewIntValue( newValues ) );
	}
}

//
// GetCurrentValue
//
// Return current values to UI
//
G4String SprayParticleGunMessenger::GetCurrentValue( G4UIcommand *command )
{
	if (command == positionCmd) {
		return positionCmd->ConvertToString( gun->GetPosition(), "cm" );
	}
	else if (command == xSprayCmd) {
		return xSprayCmd->ConvertToString( gun->GetXSpray() );
	}
	else if (command == ySprayCmd) {
		return ySprayCmd->ConvertToString( gun->GetYSpray() );
	}
	else if (command == zSprayCmd) {
		return zSprayCmd->ConvertToString( gun->GetZSpray() );
	}

	return "baloney";
}	

