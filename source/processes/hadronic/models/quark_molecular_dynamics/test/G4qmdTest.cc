#include "globals.hh"
#include "G4ios.hh"
#include "G4qmd.hh"

int main()
{

  G4String & InitialFile = "";
  G4double OutputStepping = 0.5;

  G4cout << "Please enter filename of qMD data for initialization:" << G4endl;
  G4cin >> InitialFile;

  G4cout << "-----------------------------------------" << G4endl;
  G4cout << "Please enter time steps of output:" << G4endl;
  G4cin >> OutputStepping;

	if (InitialFile == "") {
  	G4cout << "Initial file empty; using default file." << G4endl;
  }
	else {
  	G4cout << "Initial data file: " << InitialFile << G4endl;
	}



  G4qmd theRun(InitialFile);

  theRun.SetupFromFile();
  theRun.SetOutputTimestep(OutputStepping);
  theRun.justRun();


}



















