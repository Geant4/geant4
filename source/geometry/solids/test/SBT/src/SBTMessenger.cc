//
// SBTMessenger.cc
//
// Implementation of the messenger for controlling test 3
//

#include "SBTMessenger.hh"
#include "SBTrun.hh"
#include "SBTVisManager.hh"
#include "G4SolidQuery.hh"

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
SBTMessenger::SBTMessenger( const G4String prefix, const G4SolidQuery *theSolidQuery, SBTVisManager *theVisManager )
{
	//
	// Store solid query
	//
	solidQuery = theSolidQuery;
	
	//
	// Store visualization manager
	//
	visManager = theVisManager;

	//
	// Defaults (of locally stored values)
	//
	errorFile = "sbt.log";
	
	//
	// Create tester
	//
	tester = new SBTrun();
	
	//
	// Declare directory
	//
	test3Directory = new G4UIdirectory( prefix );
	test3Directory->SetGuidance( "Controls for CSG batch test" );
	
	//
	// Target command
	//
	G4String com = prefix+"target";
	targetCmd = new G4UIcmdWith3VectorAndUnit( com, this );
	targetCmd->SetGuidance( "Center of distribution of random points" );
	targetCmd->SetParameterName( "X", "Y", "Z", true, true );
	
	//
	// Widths command
	//
	com = prefix+"widths";
	widthsCmd = new G4UIcmdWith3VectorAndUnit( com, this );
	widthsCmd->SetGuidance( "Widths of distribution of random points" );
	widthsCmd->SetParameterName( "Dx", "Dy", "Dz", true, true );
	
	//
	// Grid Size command
	//
	com = prefix+"gridSizes";
	gridSizesCmd = new G4UIcmdWith3VectorAndUnit( com, this );
	gridSizesCmd->SetGuidance( "Grid size, or zero for no grid" );
	gridSizesCmd->SetParameterName( "Dx", "Dy", "Dz", true, true );
	
	//
	// Max Points command
	//
	com = prefix+"maxPoints";
	maxPointsCmd = new G4UIcmdWithAnInteger( com, this );
	maxPointsCmd->SetGuidance( "Maximum number of points before test ends" );

	//
	// Max Errors command
	//	
	com = prefix+"maxErrors";
	maxErrorsCmd = new G4UIcmdWithAnInteger( com, this );
	maxErrorsCmd->SetGuidance( "Maximum number of errors before test ends" );
	
	//
	// Error filename command
	//
	com = prefix+"errorFileName";
	errorFileCmd = new G4UIcmdWithAString( com, this );
	errorFileCmd->SetGuidance( "Filename in which to send error listings" );
	
	//
	// Run command
	//
	com = prefix+"run";
	runCmd = new G4UIcmdWithoutParameter( com, this );
	runCmd->SetGuidance( "Execute test 3" );
	
	//
	// Debug commands
	//
	com = prefix+"draw";
	drawCmd = new G4UIcmdWithAnInteger( com, this );
	drawCmd->SetGuidance( "Draw error listed in log file" );
	
	com = prefix+"debugInside";
	debugInsideCmd = new G4UIcmdWithAnInteger( com, this );
	debugInsideCmd->SetGuidance( "Call G4VSolid::Inside for error listed in log file" );
	
	com = prefix+"debugToInP";
	debugToInPCmd = new G4UIcmdWithAnInteger( com, this );
	debugToInPCmd->SetGuidance( "Call G4VSolid::DistanceToIn(p) for error listed in log file" );
	
	com = prefix+"debugToInPV";
	debugToInPVCmd = new G4UIcmdWithAnInteger( com, this );
	debugToInPVCmd->SetGuidance( "Call G4VSolid::DistanceToIn(p,v) for error listed in log file" );
	
	com = prefix+"debugToOutP";
	debugToOutPCmd = new G4UIcmdWithAnInteger( com, this );
	debugToOutPCmd->SetGuidance( "Call G4VSolid::DistanceToOut(p) for error listed in log file" );
	
	com = prefix+"debugToOutPV";
	debugToOutPVCmd = new G4UIcmdWithAnInteger( com, this );
	debugToOutPVCmd->SetGuidance( "Call G4VSolid::DistanceToOut(p,v) for error listed in log file" );

	//
	// Pause command
	//
	com = prefix+"pause";
	pauseCmd = new G4UIcmdWithoutParameter( com, this );
	pauseCmd->SetGuidance( " Wait for a key " );
	
}


//
// Destructor
//
SBTMessenger::~SBTMessenger()
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
void SBTMessenger::InvokeTest3()
{
	//
	// Is there a Solid to test?
	//
	G4VSolid *testSolid = solidQuery->GetSolid();
	if (testSolid == 0) {
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
	tester->RunTest( testSolid, logFile );
}


//
// Debug
//
// Run an event that should (approximately) duplicate the conditions
// of one of the bugs discovered in a previous test3 run and stored
// in a log file
//
void SBTMessenger::Debug( const G4int errorIndex,  SBTMessenger::Debugger *debugger )
{
	//
	// Is there a Solid to test?
	//
	// If I was really clever (I'm not) I would also check
	// to make sure we have the same solid type and even the
	// same solid parameters.
	//
	G4VSolid *testSolid = solidQuery->GetSolid();
	if (testSolid == 0) {
		G4cerr << "Please initialize geometry before debugging/drawing" << G4endl;
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
G4int SBTMessenger::DrawError::DebugMe( G4std::ifstream &logFile, const G4int errorIndex )
{
	if (visManager) {
		//
		// Make sure we're kosher
		//
		if (visManager->BuildFakeWorld()) return 2;
		
		return tester->DrawError( testSolid, logFile, errorIndex, visManager );
	}
	
	
	G4cerr << "Visualization is not available in this executable" << G4endl;
	return 1;
}

G4int SBTMessenger::DebugInside::DebugMe( G4std::ifstream &logFile, const G4int errorIndex )
{
	return tester->DebugInside( testSolid, logFile, errorIndex );
}

G4int SBTMessenger::DebugToInP::DebugMe( G4std::ifstream &logFile, const G4int errorIndex )
{
	return tester->DebugToInP( testSolid, logFile, errorIndex );
}

G4int SBTMessenger::DebugToInPV::DebugMe( G4std::ifstream &logFile, const G4int errorIndex )
{
	return tester->DebugToInPV( testSolid, logFile, errorIndex );
}

G4int SBTMessenger::DebugToOutP::DebugMe( G4std::ifstream &logFile, const G4int errorIndex )
{
	return tester->DebugToOutP( testSolid, logFile, errorIndex );
}

G4int SBTMessenger::DebugToOutPV::DebugMe( G4std::ifstream &logFile, const G4int errorIndex )
{
	return tester->DebugToOutPV( testSolid, logFile, errorIndex );
}


//
// SetNewValue
//
// Call by the UI when user requests a change
//
void SBTMessenger::SetNewValue( G4UIcommand *command, G4String newValues )
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
	else if (command == drawCmd) {
		SBTMessenger::DrawError debugger( solidQuery->GetSolid(), tester, visManager );
		Debug( drawCmd->GetNewIntValue( newValues ), &debugger );
	}
	else if (command == debugInsideCmd) {
		SBTMessenger::DebugInside debugger( solidQuery->GetSolid(), tester );
		Debug( debugInsideCmd->GetNewIntValue( newValues ), &debugger );
	}
	else if (command == debugToInPCmd) {
		SBTMessenger::DebugToInP debugger( solidQuery->GetSolid(), tester );
		Debug( debugToInPCmd->GetNewIntValue( newValues ), &debugger );
	}
	else if (command == debugToInPVCmd) {
		SBTMessenger::DebugToInPV debugger( solidQuery->GetSolid(), tester );
		Debug( debugToInPVCmd->GetNewIntValue( newValues ), &debugger );
	}
	else if (command == debugToOutPCmd) {
		SBTMessenger::DebugToOutP debugger( solidQuery->GetSolid(), tester );
		Debug( debugToOutPCmd->GetNewIntValue( newValues ), &debugger );
	}
	else if (command == debugToOutPVCmd) {
		SBTMessenger::DebugToOutPV debugger( solidQuery->GetSolid(), tester );
		Debug( debugToOutPVCmd->GetNewIntValue( newValues ), &debugger );
	}
	else if (command == pauseCmd) {
	  char c;
	  
	  G4cout << "Pause" << G4endl ;
	  G4cin.get(c);
	}
	else {
		G4Exception( "Unrecognized command" );
	}
}


//
// GetCurrentValue
//
G4String SBTMessenger::GetCurrentValue( G4UIcommand *command )
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
	else if (command == drawCmd ||
		 command == debugInsideCmd ||
		 command == debugToInPCmd ||
		 command == debugToInPVCmd ||
		 command == debugToOutPCmd ||
		 command == debugToOutPVCmd    ) {
		return "";
	}
	
	G4Exception( "Unrecognized command" );
	return "foo!";
}
	
