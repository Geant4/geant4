//
// FredVoxelTestMessenger.cc
//

#include "FredVoxelTestMessenger.hh"
#include "FredVoxelTest.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"

#include "G4VVisManager.hh"
#include "G4UImanager.hh"

//
// Constructor
//
FredVoxelTestMessenger::FredVoxelTestMessenger()
{
	voxelTest = new FredVoxelTest();
	
	//
	// Declare directory
	//
	voxelTestDirectory = new G4UIdirectory( "/fred/voxel/" );
	voxelTestDirectory->SetGuidance( "Controls for voxel test" );
	
	//
	// X min/max command
	//
	xMinMaxCmd = new G4UIcmdWith3VectorAndUnit( "/fred/voxel/xextent", this );
	xMinMaxCmd->SetGuidance( "Max/min value along voxel x axis" );
	xMinMaxCmd->SetParameterName( "Min", "Max", "NotUsed", true, true );
	
	//
	// Y min/max command
	//
	yMinMaxCmd = new G4UIcmdWith3VectorAndUnit( "/fred/voxel/yextent", this );
	yMinMaxCmd->SetGuidance( "Max/min value along voxel y axis" );
	yMinMaxCmd->SetParameterName( "Min", "Max", "NotUsed", true, true );
	
	//
	// Z min/max command
	//
	zMinMaxCmd = new G4UIcmdWith3VectorAndUnit( "/fred/voxel/zextent", this );
	zMinMaxCmd->SetGuidance( "Max/min value along voxel z axis" );
	zMinMaxCmd->SetParameterName( "Min", "Max", "NotUsed", true, true );

	//
	// Position command
	//
	posCmd = new G4UIcmdWith3VectorAndUnit( "/fred/voxel/position", this );
	posCmd->SetGuidance( "Position of voxel" );
	posCmd->SetParameterName( "X", "Y", "Z", true, true );
	
	//
	// Rotate commands
	//
	xRotateCmd = new G4UIcmdWithADouble( "/fred/voxel/xrotate", this );
	xRotateCmd->SetGuidance( "Rotate around x axis" );
	
	yRotateCmd = new G4UIcmdWithADouble( "/fred/voxel/yrotate", this );
	yRotateCmd->SetGuidance( "Rotate around y axis" );
	
	zRotateCmd = new G4UIcmdWithADouble( "/fred/voxel/zrotate", this );
	zRotateCmd->SetGuidance( "Rotate around z axis" );
	
	levelCmd = new G4UIcmdWithoutParameter( "/fred/voxel/level", this );
	levelCmd->SetGuidance( "Reset rotation" );
	
	//
	// Draw command
	//
	drawCmd = new G4UIcmdWithAString( "/fred/voxel/draw", this );
	drawCmd->SetGuidance( "Test and draw the voxel" );
	drawCmd->SetCandidates( "x y z" );
	
	//
	// Test command
	//
	testCmd = new G4UIcmdWithAString( "/fred/voxel/test", this );
	testCmd->SetGuidance( "Apply the voxel to the test volume" );
	testCmd->SetCandidates( "x y z" );
	
	//
	// Reset command
	//
	resetCmd = new G4UIcmdWithoutParameter( "/fred/voxel/reset", this );
	resetCmd->SetGuidance( "Start with a fresh voxel" );
}


//
// Destructor
//
FredVoxelTestMessenger::~FredVoxelTestMessenger()
{
	delete xMinMaxCmd;
	delete yMinMaxCmd;
	delete zMinMaxCmd;
	delete posCmd;
	delete drawCmd;
	delete testCmd;
	delete voxelTestDirectory;
	
	delete voxelTest;
}


//
// InvokeTest
//
void FredVoxelTestMessenger::InvokeTest( G4String request )
{
	if (testVolume== 0) {
		G4cerr << "Please initialize geometry before running test 3" << G4endl;
		G4cerr << "Test 3 ABORTED" << G4endl;
		return;
	}
	
	EAxis axis;
	
	if ( request == "x" ) 
		axis = kXAxis;
	else if ( request == "y" )
		axis = kYAxis;
	else if ( request == "z" )
		axis = kZAxis;
	else {
		G4cerr << "Please specify x, y, or z" << G4endl;
		return;
	}
	
	voxelTest->Test( axis, testVolume );
}


//
// Draw
//
void FredVoxelTestMessenger::Draw()
{
	if (G4VVisManager::GetConcreteInstance()) {
        	G4UImanager *UI = G4UImanager::GetUIpointer();

                // Prepare
                UI->ApplyCommand( "/vis~/clear/view" );
 
                // Draw detector
                UI->ApplyCommand( "/vis~/draw/current" );
	
		// Draw voxel
		voxelTest->Draw();
		
                // Finish
                UI->ApplyCommand( "/vis~/show/view" );
	}
}


//
// SetNewVaolume
//
// Call by the UI when user requests a change
//
void FredVoxelTestMessenger::SetNewValue( G4UIcommand *command, G4String newValues )
{
	if (command == xMinMaxCmd) {
		G4ThreeVector param = xMinMaxCmd->GetNew3VectorValue( newValues );
		voxelTest->SetExtent( kXAxis, param.x(), param.y() );
	} 
	else if (command == yMinMaxCmd) {
		G4ThreeVector param = yMinMaxCmd->GetNew3VectorValue( newValues );
		voxelTest->SetExtent( kYAxis, param.x(), param.y() );
	} 
	else if (command == zMinMaxCmd) {
		G4ThreeVector param = zMinMaxCmd->GetNew3VectorValue( newValues );
		voxelTest->SetExtent( kZAxis, param.x(), param.y() );
	} 
	else if (command == posCmd) {
		voxelTest->SetOrigin( posCmd->GetNew3VectorValue( newValues ) );
	}
	else if (command == xRotateCmd) {
		voxelTest->Rotate( kXAxis, deg*xRotateCmd->GetNewDoubleValue( newValues ) );
	}
	else if (command == yRotateCmd) {
		voxelTest->Rotate( kYAxis, deg*yRotateCmd->GetNewDoubleValue( newValues ) );
	}
	else if (command == zRotateCmd) {
		voxelTest->Rotate( kZAxis, deg*zRotateCmd->GetNewDoubleValue( newValues ) );
	}
	else if (command == levelCmd) {
		voxelTest->ResetRotation();
	}
	else if (command == drawCmd) {
		InvokeTest( newValues );
		Draw();
	}
	else if (command == testCmd) {
		InvokeTest( newValues );
	}
	else if (command == resetCmd) {
		delete voxelTest;
		voxelTest = new FredVoxelTest();
	}
	else {
		G4Exception( "Unrecognized command" );
	}
}





//
// GetCurrentValue
//
G4String FredVoxelTestMessenger::GetCurrentValue( G4UIcommand *command )
{
	return "fred";
}
	
