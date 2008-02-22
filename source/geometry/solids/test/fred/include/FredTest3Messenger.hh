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
//
// FredTest3Messenger.hh
//
// Definition of the messenger for controlling test 3
//
#ifndef FredTest3Messenger_hh
#define FredTest3Messenger_hh

#include "G4UImessenger.hh"
#include "globals.hh"
#include <fstream>

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

		virtual G4int DebugMe( std::ifstream &logFile, const G4int errorIndex ) = 0;

		protected:
		const G4VSolid *testVolume;
		const FredTest3 *tester;
	};
	
	class DebugError : public FredTest3Messenger::Debugger {
		public:
		DebugError( const G4VSolid *aVolume, const FredTest3 *aTester ) : Debugger(aVolume,aTester) {;}
		G4int DebugMe( std::ifstream &logFile, const G4int errorIndex );
	};
	class DebugInside : public FredTest3Messenger::Debugger {
		public:
		DebugInside( const G4VSolid *aVolume, const FredTest3 *aTester ) : Debugger(aVolume,aTester) {;}
		G4int DebugMe( std::ifstream &logFile, const G4int errorIndex );
	};
	class DebugToInP : public FredTest3Messenger::Debugger {
		public:
		DebugToInP( const G4VSolid *aVolume, const FredTest3 *aTester ) : Debugger(aVolume,aTester) {;}
		G4int DebugMe( std::ifstream &logFile, const G4int errorIndex );
	};
	class DebugToInPV : public FredTest3Messenger::Debugger {
		public:
		DebugToInPV( const G4VSolid *aVolume, const FredTest3 *aTester ) : Debugger(aVolume,aTester) {;}
		G4int DebugMe( std::ifstream &logFile, const G4int errorIndex );
	};
	class DebugToOutP : public FredTest3Messenger::Debugger {
		public:
		DebugToOutP( const G4VSolid *aVolume, const FredTest3 *aTester ) : Debugger(aVolume,aTester) {;}
		G4int DebugMe( std::ifstream &logFile, const G4int errorIndex );
	};
	class DebugToOutPV : public FredTest3Messenger::Debugger  {
		public:
		DebugToOutPV( const G4VSolid *aVolume, const FredTest3 *aTester ) : Debugger(aVolume,aTester) {;}
		G4int DebugMe( std::ifstream &logFile, const G4int errorIndex );
	};
	
        void InvokeTest3();
        void InvokeRunDebug();
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
	G4UIcmdWithoutParameter		*runDebugCmd;
	G4UIcmdWithAnInteger		*debugCmd;
	G4UIcmdWithAnInteger		*debugInsideCmd;
	G4UIcmdWithAnInteger		*debugToInPCmd;
	G4UIcmdWithAnInteger		*debugToInPVCmd;
	G4UIcmdWithAnInteger		*debugToOutPCmd;
	G4UIcmdWithAnInteger		*debugToOutPVCmd;
};

#endif
