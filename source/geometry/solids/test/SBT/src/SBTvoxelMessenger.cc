//
// SBTvoxelMessenger.cc
//
// Implementation of the messenger for the voxel test
//

#include "SBTvoxelMessenger.hh"
#include "SBTvoxel.hh"
#include "SBTVisManager.hh"
#include "G4SolidQuery.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "G4ios.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3Vector.hh"

#include "G4UIcmdWithPargs.hh"
#include "G4UIcmdPargDouble.hh"
#include "G4UIcmdPargInteger.hh"

#include <fstream.h>

//
// Constructor
//
SBTvoxelMessenger::SBTvoxelMessenger( const G4String prefix, 
				      const G4SolidQuery *theSolidQuery, SBTVisManager *theVisManager )
			: tester(), errorFile("sbtvoxel.log")
{
	//
	// Initialize
	//
	visManager = theVisManager;
	translate = new G4ThreeVector;
	rotate = new G4RotationMatrix;
	voxel = 0;
	point = 0;
	limits[0] = +1000.0;
	limits[1] = -1000.0;
	axis = kXAxis;
	
	//
	// Store solid query
	//
	solidQuery = theSolidQuery;

	//
	// Declare directory
	//
	voxelDirectory = new G4UIdirectory( prefix );
	voxelDirectory->SetGuidance( "Controls for CSG batch voxel tester" );
	
	//
	// Target command
	//
	G4String command = prefix+"target";
	targetCmd = new G4UIcmdWith3VectorAndUnit( command, this );
	targetCmd->SetGuidance( "Center of distribution of random points" );
	targetCmd->SetParameterName( "X", "Y", "Z", true, true );
	
	//
	// Widths command
	//
	command = prefix+"widths";
	widthsCmd = new G4UIcmdWith3VectorAndUnit( command, this );
	widthsCmd->SetGuidance( "Widths of distribution of random points" );
	widthsCmd->SetParameterName( "Dx", "Dy", "Dz", true, true );

	//
	// Max Voxels command
	//
	command = prefix+"maxVoxels";
	maxVoxelsCmd = new G4UIcmdWithAnInteger( command, this );
	maxVoxelsCmd->SetGuidance( "Maximum number of Voxels before test ends" );

	//
	// Max Errors command
	//	
	command = prefix+"maxErrors";
	maxErrorsCmd = new G4UIcmdWithAnInteger( command, this );
	maxErrorsCmd->SetGuidance( "Maximum number of errors before test ends" );
	
	//
	// Error filename command
	//
	command = prefix+"errorFileName";
	errorFileCmd = new G4UIcmdWithAString( command, this );
	errorFileCmd->SetGuidance( "Filename in which to send error listings" );
	
	//
	// Run command
	//
	command = prefix+"run";
	runCmd = new G4UIcmdWithoutParameter( command, this );
	runCmd->SetGuidance( "Execute voxel test" );
	
	//
	// Picture directory
	//
	command = prefix+"picture/";
	pictDirectory = new G4UIdirectory( command );
	pictDirectory->SetGuidance( "Controls for drawing CSG voxel tests" );
	
	//
	// picture/voxel command
	//
	pictVoxelArgs[0] = new G4UIcmdPargDouble( "xmin", 1.0, mm );
	pictVoxelArgs[1] = new G4UIcmdPargDouble( "xmax", 1.0, mm );
	pictVoxelArgs[2] = new G4UIcmdPargDouble( "ymin", 1.0, mm );
	pictVoxelArgs[3] = new G4UIcmdPargDouble( "ymax", 1.0, mm );
	pictVoxelArgs[4] = new G4UIcmdPargDouble( "zmin", 1.0, mm );
	pictVoxelArgs[5] = new G4UIcmdPargDouble( "zmax", 1.0, mm );
	command = prefix+"picture/voxel";
	pictVoxelCmd = new G4UIcmdWithPargs( command, this, pictVoxelArgs, 6 );
	pictVoxelCmd->SetGuidance( "Specifies the voxel to draw" );

	//
	// picture/translate command
	//
	command = prefix+"picture/translate";
	pictTranCmd = new G4UIcmdWith3Vector( command, this );
	pictTranCmd->SetGuidance( "Specifies the voxel translation" );
	
	//
	// picture/rotate command
	//
	pictRotArgs[0] = new G4UIcmdPargDouble( "xaxis", 1.0, 1 );
	pictRotArgs[1] = new G4UIcmdPargDouble( "yaxis", 1.0, 1 );
	pictRotArgs[2] = new G4UIcmdPargDouble( "zaxis", 1.0, 1 );
	pictRotArgs[3] = new G4UIcmdPargDouble( "amount", 1.0, 1 );
	command = prefix+"picture/rotate";
	pictRotCmd = new G4UIcmdWithPargs( command, this, pictRotArgs, 4 );
	pictRotCmd->SetGuidance( "Specifies rotation matrix" );

	//
	// picture/point command
	//
	command = prefix+"picture/point";
	pictPointCmd = new G4UIcmdWith3Vector( command, this );
	pictPointCmd->SetGuidance( "Specifies one point to include in the picture" );
	
	//
	// picture/limit command
	//
	pictLimitArgs[0] = new G4UIcmdPargInteger( "axis", 1 );
	pictLimitArgs[1] = new G4UIcmdPargDouble( "min", 1.0, mm );
	pictLimitArgs[2] = new G4UIcmdPargDouble( "max", 1.0, mm );
	command = prefix+"picture/limit";
	pictLimitCmd = new G4UIcmdWithPargs( command, this, pictLimitArgs, 3 );
	pictLimitCmd->SetGuidance( "Specifies voxel limits" );
	
	//
	// picture/draw command
	//
	command = prefix+"picture/draw";
	pictDrawCmd = new G4UIcmdWithoutParameter( command, this );
	pictDrawCmd->SetGuidance( "Draw the picture" );

	//
	// picture/debug command
	//
	command = prefix+"picture/debug";
	pictDebugCmd = new G4UIcmdWithoutParameter( command, this );
	pictDebugCmd->SetGuidance( "Call CalculateExtent for the picture" );
}



//
// Destructor
//
SBTvoxelMessenger::~SBTvoxelMessenger()
{
	if (voxel) delete voxel;
	if (point) delete point;

//	delete targetCmd;
//	delete widthsCmd;
//	delete maxVoxelsCmd;
//	delete maxVoxelsCmd;
//	delete errorFileCmd;
	delete voxelDirectory;
	
	delete pictDirectory;
}
	

//
// InvokeTest
//
// Run the test
//
void SBTvoxelMessenger::InvokeTest()
{
	//
	// Is there a Solid to test?
	//
	G4VSolid *testSolid = solidQuery->GetSolid();
	if (testSolid == 0) {
		G4cerr << "Please initialize geometry before running test" << endl;
		G4cerr << "Test ABORTED" << endl;
		return;
	}

	//
	// Open output file
	//
	ofstream logFile( errorFile );
	
	//
	// Run the test
	//
	tester.RunTest( testSolid, logFile );
}


void SBTvoxelMessenger::Draw()
{
	if (!visManager) 
		G4cerr << "Visualization is not available in this executable" << endl;

	//
	// Is there a Solid to test?
	//
	G4VSolid *testSolid = solidQuery->GetSolid();
	if (testSolid == 0) {
		G4cerr << "Please initialize geometry before running test" << endl;
		return;
	}
	
	//
	// Is there a voxel to test?
	//
	if (voxel == 0) {
		G4cerr << "Please specify a voxel to draw using /pict/voxel command" << endl;
		return;
	}

	//
	// Go do it
	//
	visManager->BuildFakeWorld();
	
	G4AffineTransform transform( *rotate, *translate );
	
	tester.Draw( testSolid, *voxel, transform, point, point ? 1 : 0, axis, limits, visManager );
}


void SBTvoxelMessenger::Debug()
{
	//
	// Is there a Solid to test?
	//
	G4VSolid *testSolid = solidQuery->GetSolid();
	if (testSolid == 0) {
		G4cerr << "Please initialize geometry before running test" << endl;
		return;
	}
	
	//
	// Is there a voxel to test?
	//
	if (voxel == 0) {
		G4cerr << "Please specify a voxel to debug using /pict/voxel command" << endl;
		return;
	}

	//
	// Go do it
	//
	G4AffineTransform transform( *rotate, *translate );
	
	tester.Debug( testSolid, axis, *voxel, transform, point );
}


//
// SetNewValue
//
// Call by the UI when user requests a change
//
void SBTvoxelMessenger::SetNewValue( G4UIcommand *command, G4String newValues )
{
	if (command == targetCmd) {
		tester.SetTarget( targetCmd->GetNew3VectorValue( newValues ) );
	}
	else if (command == widthsCmd) {
		tester.SetWidths( widthsCmd->GetNew3VectorValue( newValues ) );
	}
	else if (command == maxVoxelsCmd) {
		tester.SetMaxVoxels( maxVoxelsCmd->GetNewIntValue( newValues ) );
	}
	else if (command == maxErrorsCmd) {
		tester.SetMaxErrors( maxErrorsCmd->GetNewIntValue( newValues ) );
	}
	else if (command == errorFileCmd) {
		errorFile = newValues;
	}
	else if (command == runCmd) {
		InvokeTest();
	}
	else if (command == pictVoxelCmd) {
		if (pictVoxelCmd->GetArguments( newValues )) {
			G4UIcmdPargDouble *minX = (G4UIcmdPargDouble *)pictVoxelArgs[0],
					  *maxX = (G4UIcmdPargDouble *)pictVoxelArgs[1],
					  *minY = (G4UIcmdPargDouble *)pictVoxelArgs[2],
					  *maxY = (G4UIcmdPargDouble *)pictVoxelArgs[3],
					  *minZ = (G4UIcmdPargDouble *)pictVoxelArgs[4],
					  *maxZ = (G4UIcmdPargDouble *)pictVoxelArgs[5];

			if (voxel) delete voxel;

			voxel = new G4VoxelLimits;

			if (minX->GetValue() < maxX->GetValue())
				voxel->AddLimit( kXAxis, minX->GetValue(), maxX->GetValue() );
			if (minY->GetValue() < maxY->GetValue())
				voxel->AddLimit( kYAxis, minY->GetValue(), maxY->GetValue() );
			if (minZ->GetValue() < maxZ->GetValue())
				voxel->AddLimit( kZAxis, minZ->GetValue(), maxZ->GetValue() );
		}
		else
			G4cerr << "Syntax error" << endl;
	}
	else if (command == pictTranCmd) {
		*translate = pictTranCmd->GetNew3VectorValue( newValues );
	}
	else if (command == pictRotCmd) {
		if (pictRotCmd->GetArguments( newValues )) {
			G4UIcmdPargDouble *rotX = (G4UIcmdPargDouble *)pictRotArgs[0],
					  *rotY = (G4UIcmdPargDouble *)pictRotArgs[1],
					  *rotZ = (G4UIcmdPargDouble *)pictRotArgs[2],
					  *amount = (G4UIcmdPargDouble *)pictRotArgs[3];

			G4ThreeVector rotAxis( rotX->GetValue(), rotY->GetValue(), rotZ->GetValue() );

			delete rotate;
			rotate = new G4RotationMatrix();
			rotate->rotate( amount->GetValue(), &rotAxis );
		}
		else
			G4cerr << "Syntax error" << endl;
	}
	else if (command == pictPointCmd) {
		if (point) delete point;
		point = new G4ThreeVector();
		*point = pictPointCmd->GetNew3VectorValue( newValues );
	}
	else if (command == pictLimitCmd) {
		if (pictLimitCmd->GetArguments( newValues )) {
			static const EAxis axes[3] = { kXAxis, kYAxis, kZAxis };
			
			G4UIcmdPargInteger *axisArg = (G4UIcmdPargInteger *)pictLimitArgs[0];
			G4UIcmdPargDouble     *min  = (G4UIcmdPargDouble  *)pictLimitArgs[1],
					      *max  = (G4UIcmdPargDouble  *)pictLimitArgs[2];
			
			G4int axisIndex = axisArg->GetValue();	  
			if (axisIndex < 0 || axisIndex > 2) {
				G4cerr << "Axis value must be from 0 to 2" << endl;
			}
			else {
				axis = axes[axisIndex];
				limits[0] = min->GetValue();
				limits[1] = max->GetValue();
			}
		}
		else
			G4cerr << "Syntax error" << endl;
	}
	else if (command == pictDrawCmd) {
		Draw();
	}
	else if (command == pictDebugCmd) {
		Debug();
	}
	else {
		G4Exception( "Unrecognized command" );
	}
}

	

//
// GetCurrentValue
//
G4String SBTvoxelMessenger::GetCurrentValue( G4UIcommand *command )
{
	if (command == targetCmd) {
		return targetCmd->ConvertToString( tester.GetTarget(), "m" );
	}
	else if (command == widthsCmd) {
		return widthsCmd->ConvertToString( tester.GetWidths(), "m" );
	}
	else if (command == maxVoxelsCmd) {
		return maxVoxelsCmd->ConvertToString( tester.GetMaxVoxels() );
	}
	else if (command == maxErrorsCmd) {
		return maxErrorsCmd->ConvertToString( tester.GetMaxErrors() );
	}
	else if (command == errorFileCmd) {
		return errorFile;
	}
	else if (command == runCmd) {
		return "";
	}
	else if (command == pictVoxelCmd) {
		return "";
	}
	else if (command == pictTranCmd) {
		return "";
	}
	else if (command == pictRotCmd) {
		return "";
	}
	else if (command == pictPointCmd) {
		return "";
	}
	else if (command == pictLimitCmd) {
		return "";
	}
	else if (command == pictDebugCmd) {
		return "";
	}
	else if (command == pictDrawCmd) {
		return "";
	}
	
	G4Exception( "Unrecognized command" );
	return "foo!";
}
