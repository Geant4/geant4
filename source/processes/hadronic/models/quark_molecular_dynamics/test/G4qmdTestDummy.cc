#include "globals.hh"
#include "G4ios.hh"
#include "G4qmdDummy.hh"

int main()
{

  G4String & InitialFile = "";

  G4cout << "Please enter filename of qMD data for initialization:" << G4endl;
  G4cin >> InitialFile;

	if (InitialFile == "") {
  	G4cout << "Initial file empty; using default file." << G4endl;
  }
	else {
  	G4cout << "Initial data file: " << InitialFile << G4endl;
	}

  G4qmdDummy theRun(InitialFile);
	
//  theRun.SetupFromFile();
//  theRun.justRun();
  G4cout << "   ...using inital file " << theRun.GetInputFile() << G4endl;

}



















