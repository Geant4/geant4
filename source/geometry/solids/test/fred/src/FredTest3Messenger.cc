//
// FredTest3Messenger.cc
//
// Implementation of the messenger for controlling test 3
//

#include "FredTest3Messenger.hh"
#include "FredTest3.hh"

#include "G4ios.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"

#include "g4std/fstream"

//
// Constructor
//
FredTest3Messenger::FredTest3Messenger()
{
	//
	// Defaults (of locally stored values)
	//
	errorFile = "test3.log";
	
	//
	// Create tester
	//
	tester = new FredTest3();
	
	//
	// Declare directory
	//
	test3Directory = new G4UIdirectory( "/fred/test3/" );
	test3Directory->SetGuidance( "Controls for volume test 3" );
	
	//
	// Target command
	//
	targetCmd = new G4UIcmdWith3VectorAndUnit( "/fred/test3/target", this );
	targetCmd->SetGuidance( "Center of distribution of random points" );
	targetCmd->SetParameterName( "X", "Y", "Z", true, true );
	
	//
	// Widths command
	//
	widthsCmd = new G4UIcmdWith3VectorAndUnit( "/fred/test3/widths", this );
	widthsCmd->SetGuidance( "Widths of distribution of random points" );
	widthsCmd->SetParameterName( "Dx", "Dy", "Dz", true, true );
	
	//
	// Grid Size command
	//
	gridSizesCmd = new G4UIcmdWith3VectorAndUnit( "/fred/test3/gridSizes", this );
	gridSizesCmd->SetGuidance( "Grid size, or zero for no grid" );
	gridSizesCmd->SetParameterName( "Dx", "Dy", "Dz", true, true );
	
	//
	// Max Points command
	//
	maxPointsCmd = new G4UIcmdWithAnInteger( "/fred/test3/maxPoints", this );
	maxPointsCmd->SetGuidance( "Maximum number of points before test ends" );

	//
	// Max Errors command
	//	
	maxErrorsCmd = new G4UIcmdWithAnInteger( "/fred/test3/maxErrors", this );
	maxErrorsCmd->SetGuidance( "Maximum number of errors before test ends" );
	
	//
	// Error filename command
	//
	errorFileCmd = new G4UIcmdWithAString( "/fred/test3/errorFileName", this );
	errorFileCmd->SetGuidance( "Filename in which to send error listings" );
	
	//
	// Run command
	//
	runCmd = new G4UIcmdWithoutParameter( "/fred/test3/run", this );
	runCmd->SetGuidance( "Execute test 3" );
	
	//
	// Debug commands
	//
	debugCmd = new G4UIcmdWithAnInteger( "/fred/test3/debug", this );
	debugCmd->SetGuidance( "Debug error listed in log file" );
	
	debugInsideCmd = new G4UIcmdWithAnInteger( "/fred/test3/debugInside", this );
	debugInsideCmd->SetGuidance( "Call G4VSolid::Inside for error listed in log file" );
	
	debugToInPCmd = new G4UIcmdWithAnInteger( "/fred/test3/debugToInP", this );
	debugToInPCmd->SetGuidance( "Call G4VSolid::DistanceToIn(p) for error listed in log file" );
	
	debugToInPVCmd = new G4UIcmdWithAnInteger( "/fred/test3/debugToInPV", this );
	debugToInPVCmd->SetGuidance( "Call G4VSolid::DistanceToIn(p,v) for error listed in log file" );
	
	debugToOutPCmd = new G4UIcmdWithAnInteger( "/fred/test3/debugToOutP", this );
	debugToOutPCmd->SetGuidance( "Call G4VSolid::DistanceToOut(p) for error listed in log file" );
	
	debugToOutPVCmd = new G4UIcmdWithAnInteger( "/fred/test3/debugToOutPV", this );
	debugToOutPVCmd->SetGuidance( "Call G4VSolid::DistanceToOut(p,v) for error listed in log file" );
}


//
// Destructor
//
FredTest3Messenger::~FredTest3Messenger()
{
	delete targetCmd;
	delete widthsCmd;
	delete gridSizesCmd;
	delete maxPointsCmd;
	delete maxErrorsCmd;
	delete errorFileCmd;
	delete test3Directory;
}


//
// InvokeTest3
//
// Run test 3
//
void FredTest3Messenger::InvokeTest3()
{
	//
	// Is there a volume to test?
	//
	if (testVolume == 0) {
		G4cerr << "Please initialize geometry before running test 3" << G4endl;
		G4cerr << "Test 3 ABORTED" << G4endl;
		return;
	}

	//
	// Open output file
	//
	G4std::ofstream logFile( errorFile );
	
	//
	// Run the test
	//
	tester->RunTest( testVolume, logFile );
}


//
// DebugError
//
// Run an event that should (approximately) duplicate the conditions
// of one of the bugs discovered in a previous test3 run and stored
// in a log file
//
void FredTest3Messenger::Debug( const G4int errorIndex,  FredTest3Messenger::Debugger *debugger )
{
	//
	// Is there a volume to test?
	//
	// If I was really clever (I'm not) I would also check
	// to make sure we have the same volume type and even the
	// same volume parameters.
	//
	if (testVolume == 0) {
		G4cerr << "Please initialize geometry before running test 3" << G4endl;
		return;
	}

	//
	// Open output file
	//
	G4std::ifstream logFile( errorFile );
	if (!logFile) {
		G4cerr << "Cannot open input file " << errorFile << G4endl;
		return;
	}
	
	//
	// Debug
	//
	if (debugger->DebugMe( logFile, errorIndex )) {
		G4cerr << "Error reading index " << errorIndex << " from input file " << errorFile << G4endl;
	}
}


//
// DebugMe (various classes)
//
G4int FredTest3Messenger::DebugError::DebugMe( G4std::ifstream &logFile, const G4int errorIndex )
{
	return tester->DebugError( testVolume, logFile, errorIndex );
}

G4int FredTest3Messenger::DebugInside::DebugMe( G4std::ifstream &logFile, const G4int errorIndex )
{
	return tester->DebugInside( testVolume, logFile, errorIndex );
}

G4int FredTest3Messenger::DebugToInP::DebugMe( G4std::ifstream &logFile, const G4int errorIndex )
{
	return tester->DebugToInP( testVolume, logFile, errorIndex );
}

G4int FredTest3Messenger::DebugToInPV::DebugMe( G4std::ifstream &logFile, const G4int errorIndex )
{
	return tester->DebugToInPV( testVolume, logFile, errorIndex );
}

G4int FredTest3Messenger::DebugToOutP::DebugMe( G4std::ifstream &logFile, const G4int errorIndex )
{
	return tester->DebugToOutP( testVolume, logFile, errorIndex );
}

G4int FredTest3Messenger::DebugToOutPV::DebugMe( G4std::ifstream &logFile, const G4int errorIndex )
{
	return tester->DebugToOutPV( testVolume, logFile, errorIndex );
}


//
// SetNewValue
//
// Call by the UI when user requests a change
//
void FredTest3Messenger::SetNewValue( G4UIcommand *command, G4String newValues )
{
	if (command == targetCmd) {
		tester->SetTarget( targetCmd->GetNew3VectorValue( newValues ) );
	}
	else if (command == widthsCmd) {
		tester->SetWidths( widthsCmd->GetNew3VectorValue( newValues ) );
	}
	else if (command == gridSizesCmd) {
		tester->SetGrids( gridSizesCmd->GetNew3VectorValue( newValues ) );
	}
	else if (command == maxPointsCmd) {
		tester->SetMaxPoints( maxPointsCmd->GetNewIntValue( newValues ) );
	}
	else if (command == maxErrorsCmd) {
		tester->SetMaxErrors( maxErrorsCmd->GetNewIntValue( newValues ) );
	}
	else if (command == errorFileCmd) {
		errorFile = newValues;
	}
	else if (command == runCmd) {
		InvokeTest3();
	}
	else if (command == debugCmd) {
		FredTest3Messenger::DebugError debugger( testVolume, tester );
		Debug( debugCmd->GetNewIntValue( newValues ), &debugger );
	}
	else if (command == debugInsideCmd) {
		FredTest3Messenger::DebugInside debugger( testVolume, tester );
		Debug( debugInsideCmd->GetNewIntValue( newValues ), &debugger );
	}
	else if (command == debugToInPCmd) {
		FredTest3Messenger::DebugToInP debugger( testVolume, tester );
		Debug( debugToInPCmd->GetNewIntValue( newValues ), &debugger );
	}
	else if (command == debugToInPVCmd) {
		FredTest3Messenger::DebugToInPV debugger( testVolume, tester );
		Debug( debugToInPVCmd->GetNewIntValue( newValues ), &debugger );
	}
	else if (command == debugToOutPCmd) {
		FredTest3Messenger::DebugToOutP debugger( testVolume, tester );
		Debug( debugToOutPCmd->GetNewIntValue( newValues ), &debugger );
	}
	else if (command == debugToOutPVCmd) {
		FredTest3Messenger::DebugToOutPV debugger( testVolume, tester );
		Debug( debugToOutPVCmd->GetNewIntValue( newValues ), &debugger );
	}
	else {
		G4Exception( "Unrecognized command" );
	}
}


//
// GetCurrentValue
//
G4String FredTest3Messenger::GetCurrentValue( G4UIcommand *command )
{
	if (command == targetCmd) {
		return targetCmd->ConvertToString( tester->GetTarget(), "m" );
	}
	else if (command == widthsCmd) {
		return widthsCmd->ConvertToString( tester->GetWidths(), "m" );
	}
	else if (command == gridSizesCmd) {
		return gridSizesCmd->ConvertToString( tester->GetGrids(), "m" );
	}
	else if (command == maxPointsCmd) {
		return maxPointsCmd->ConvertToString( tester->GetMaxPoints() );
	}
	else if (command == maxErrorsCmd) {
		return maxErrorsCmd->ConvertToString( tester->GetMaxErrors() );
	}
	else if (command == errorFileCmd) {
		return errorFile;
	}
	else if (command == runCmd) {
		return "";
	}
	else if (command == debugCmd ||
		 command == debugInsideCmd ||
		 command == debugToInPCmd ||
		 command == debugToInPVCmd ||
		 command == debugToOutPCmd ||
		 command == debugToOutPVCmd    ) {
		return "";
	}
	
	G4Exception( "Unrecognized command" );
	return "";
}
	
