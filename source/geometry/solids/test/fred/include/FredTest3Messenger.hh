//
// FredTest3Messenger.hh
//
// Definition of the messenger for controlling test 3
//
#ifndef FredTest3Messenger_hh
#define FredTest3Messenger_hh

#include "G4UImessenger.hh"
#include "globals.hh"
#include "g4std/fstream"

class FredTest3;
class G4VSolid;

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

class FredTest3Messenger : public G4UImessenger
{
	public:
	FredTest3Messenger();
	~FredTest3Messenger();
	
	void SetNewValue( G4UIcommand *command, G4String newValues );
	G4String GetCurrentValue( G4UIcommand *command );
	
	inline void SetTestVolume( const G4VSolid *newTestVolume ) { testVolume = newTestVolume; }
	inline const G4VSolid *GetTestVolume() const { return testVolume; }

	protected:
	class Debugger {
		public:
		Debugger( const G4VSolid *aVolume, const FredTest3 *aTester ) 
					{ testVolume = aVolume; tester = aTester; }
		virtual ~Debugger() {;}

		virtual G4int DebugMe( G4std::ifstream &logFile, const G4int errorIndex ) = 0;

		protected:
		const G4VSolid *testVolume;
		const FredTest3 *tester;
	};
	
	class DebugError : public FredTest3Messenger::Debugger {
		public:
		DebugError( const G4VSolid *aVolume, const FredTest3 *aTester ) : Debugger(aVolume,aTester) {;}
		G4int DebugMe( G4std::ifstream &logFile, const G4int errorIndex );
	};
	class DebugInside : public FredTest3Messenger::Debugger {
		public:
		DebugInside( const G4VSolid *aVolume, const FredTest3 *aTester ) : Debugger(aVolume,aTester) {;}
		G4int DebugMe( G4std::ifstream &logFile, const G4int errorIndex );
	};
	class DebugToInP : public FredTest3Messenger::Debugger {
		public:
		DebugToInP( const G4VSolid *aVolume, const FredTest3 *aTester ) : Debugger(aVolume,aTester) {;}
		G4int DebugMe( G4std::ifstream &logFile, const G4int errorIndex );
	};
	class DebugToInPV : public FredTest3Messenger::Debugger {
		public:
		DebugToInPV( const G4VSolid *aVolume, const FredTest3 *aTester ) : Debugger(aVolume,aTester) {;}
		G4int DebugMe( G4std::ifstream &logFile, const G4int errorIndex );
	};
	class DebugToOutP : public FredTest3Messenger::Debugger {
		public:
		DebugToOutP( const G4VSolid *aVolume, const FredTest3 *aTester ) : Debugger(aVolume,aTester) {;}
		G4int DebugMe( G4std::ifstream &logFile, const G4int errorIndex );
	};
	class DebugToOutPV : public FredTest3Messenger::Debugger  {
		public:
		DebugToOutPV( const G4VSolid *aVolume, const FredTest3 *aTester ) : Debugger(aVolume,aTester) {;}
		G4int DebugMe( G4std::ifstream &logFile, const G4int errorIndex );
	};
	
        void InvokeTest3();
        void Debug( const G4int errorIndex, FredTest3Messenger::Debugger *debugger );
	
	private:
	FredTest3	*tester;
	const G4VSolid	*testVolume;
	
	G4String	errorFile;
	
	G4UIdirectory			*test3Directory;
	G4UIcmdWith3VectorAndUnit	*targetCmd;
	G4UIcmdWith3VectorAndUnit	*widthsCmd;
	G4UIcmdWith3VectorAndUnit	*gridSizesCmd;
	G4UIcmdWithAnInteger		*maxPointsCmd;
	G4UIcmdWithAnInteger		*maxErrorsCmd;
	G4UIcmdWithAString		*errorFileCmd;
	G4UIcmdWithoutParameter		*runCmd;
	G4UIcmdWithAnInteger		*debugCmd;
	G4UIcmdWithAnInteger		*debugInsideCmd;
	G4UIcmdWithAnInteger		*debugToInPCmd;
	G4UIcmdWithAnInteger		*debugToInPVCmd;
	G4UIcmdWithAnInteger		*debugToOutPCmd;
	G4UIcmdWithAnInteger		*debugToOutPVCmd;
};

#endif
