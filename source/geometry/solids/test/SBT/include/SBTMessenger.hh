//
// SBTMessenger.hh
//
// Definition of the messenger for controlling test 3
//
#ifndef SBTMessenger_hh
#define SBTMessenger_hh

#include "G4UImessenger.hh"
#include "globals.hh"
#include "g4std/fstream"

class SBTrun;
class SBTVisManager;

class G4VSolid;
class G4SolidQuery;

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

class SBTMessenger : public G4UImessenger
{
	public:
	SBTMessenger( const G4String prefix, const G4SolidQuery *solidQuery, SBTVisManager *visManager );
	~SBTMessenger();
	
	void SetNewValue( G4UIcommand *command, G4String newValues );
	G4String GetCurrentValue( G4UIcommand *command );
	
	inline const G4SolidQuery *GetSolidQuery() const { return solidQuery; }
	
	inline SBTVisManager *GetVisManager() const { return visManager; }
	inline void SetVisManager( SBTVisManager *theVisManager ) { visManager = theVisManager; }

	public:
	class Debugger {
		public:
		Debugger( const G4VSolid *aSolid, const SBTrun *aTester ) 
					{ testSolid = aSolid; tester = aTester; }
		virtual ~Debugger() {;}

		virtual G4int DebugMe( G4std::ifstream &logFile, const G4int errorIndex ) = 0;

		protected:
		const G4VSolid *testSolid;
		const SBTrun *tester;
	};
	
	protected:
	class DrawError : public SBTMessenger::Debugger {
		public:
		DrawError( const G4VSolid *aSolid, const SBTrun *aTester, SBTVisManager *theVisManager )
				 : Debugger(aSolid,aTester) { visManager = theVisManager;}
		G4int DebugMe( G4std::ifstream &logFile, const G4int errorIndex );
		
		protected:
		SBTVisManager *visManager;
	};
	class DebugInside : public SBTMessenger::Debugger {
		public:
		DebugInside( const G4VSolid *aSolid, const SBTrun *aTester ) : Debugger(aSolid,aTester) {;}
		G4int DebugMe( G4std::ifstream &logFile, const G4int errorIndex );
	};
	class DebugToInP : public SBTMessenger::Debugger {
		public:
		DebugToInP( const G4VSolid *aSolid, const SBTrun *aTester ) : Debugger(aSolid,aTester) {;}
		G4int DebugMe( G4std::ifstream &logFile, const G4int errorIndex );
	};
	class DebugToInPV : public SBTMessenger::Debugger {
		public:
		DebugToInPV( const G4VSolid *aSolid, const SBTrun *aTester ) : Debugger(aSolid,aTester) {;}
		G4int DebugMe( G4std::ifstream &logFile, const G4int errorIndex );
	};
	class DebugToOutP : public SBTMessenger::Debugger {
		public:
		DebugToOutP( const G4VSolid *aSolid, const SBTrun *aTester ) : Debugger(aSolid,aTester) {;}
		G4int DebugMe( G4std::ifstream &logFile, const G4int errorIndex );
	};
	class DebugToOutPV : public SBTMessenger::Debugger  {
		public:
		DebugToOutPV( const G4VSolid *aSolid, const SBTrun *aTester ) : Debugger(aSolid,aTester) {;}
		G4int DebugMe( G4std::ifstream &logFile, const G4int errorIndex );
	};
	
        void InvokeTest3();
        void Debug( const G4int errorIndex, SBTMessenger::Debugger *debugger );
	
	private:
	SBTrun		*tester;
	const G4SolidQuery *solidQuery;
	SBTVisManager	*visManager;
	
	G4String	errorFile;
	
	G4UIdirectory			*test3Directory;
	G4UIcmdWith3VectorAndUnit	*targetCmd;
	G4UIcmdWith3VectorAndUnit	*widthsCmd;
	G4UIcmdWith3VectorAndUnit	*gridSizesCmd;
	G4UIcmdWithAnInteger		*maxPointsCmd;
	G4UIcmdWithAnInteger		*maxErrorsCmd;
	G4UIcmdWithAString		*errorFileCmd;
	G4UIcmdWithoutParameter		*runCmd;
	G4UIcmdWithAnInteger		*drawCmd;
	G4UIcmdWithAnInteger		*debugInsideCmd;
	G4UIcmdWithAnInteger		*debugToInPCmd;
	G4UIcmdWithAnInteger		*debugToInPVCmd;
	G4UIcmdWithAnInteger		*debugToOutPCmd;
	G4UIcmdWithAnInteger		*debugToOutPVCmd;
	G4UIcmdWithoutParameter		*pauseCmd;
};

#endif
