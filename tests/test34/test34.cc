using namespace std;


//std includes 
#include <algorithm>
#include <iostream>
#include "G4Timer.hh"
// G4 includes 
#include "G4ios.hh"
#include "G4GlobalFastSimulationManager.hh"
#include "G4Timer.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

// my project 
#include "Tst34DetectorConstruction.hh"
#include "Tst34PhysicsList.hh"
#include "Tst34PrimaryGeneratorAction.hh"
#include "Tst34EventAction.hh"
#include "Tst34RunAction.hh"

#ifdef G4_SOLVE_TEMPLATES
#ifdef G4VIS_USE
#define G4_SOLVE_VIS_TEMPLATES
#endif
#endif
G4Timer Timer;
G4Timer Timerintern;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{ 	
	// Timer to see GFlash performance
	Timer.Start();
	
	cout<<"+-------------------------------------------------------+"<<endl;
	cout<<"|                                                       |"<<endl;
	cout<<"|          This is an example of Shower                 |"<<endl;
	cout<<"|          Parameterization with GFLASH                 |"<<endl;
	cout<<"+-------------------------------------------------------+"<<endl;
	
	G4RunManager* runManager = new G4RunManager;
	
	// UserInitialization classes (mandatory)
	cout<<"# GFlash Example: Detector Construction"<<endl;    
	runManager->SetUserInitialization(new Tst34DetectorConstruction);
	cout<<"# GFlash Example: Physics list"<<endl;
	runManager->SetUserInitialization(new Tst34PhysicsList);
	cout<<"# GFlash Example: Primary Generator"<<endl;
	runManager->SetUserAction(new Tst34PrimaryGeneratorAction);
	cout<<"# GFlash Example: User Action Classes"<<endl;
	runManager->SetUserAction(new Tst34EventAction);
	runManager->SetUserAction(new Tst34RunAction);
	
	
	
	G4UImanager* UI = G4UImanager::GetUIpointer();
	UI->ApplyCommand("/run/verbose 0");
	runManager->Initialize();
	UI->ApplyCommand("/Step/Verbose 0");
	
	if (argc==1)   // Define UI terminal for interactive mode  
	{ 
		G4UIsession * session = new G4UIterminal(new G4UItcsh);	
		//	UI->ApplyCommand("/control/execute test.mac");    
		session->SessionStart();
		//delete session;
	}
	else           // Batch mode
	{ 
		G4String s=*(argv+1);
		UI->ApplyCommand("/control/execute "+s);
	}
	
	
	delete runManager;
	
	Timer.Stop();
	cout << endl;
	cout << "******************************************";
	cout << endl;
	cout << "Total Real Elapsed Time is: "<< Timer.GetRealElapsed();
	cout << endl;
	cout << "Total System Elapsed Time: " << Timer.GetSystemElapsed();
	cout << endl;
	cout << "Total GetUserElapsed Time: " << Timer.GetUserElapsed();
	cout << endl;
	cout << "******************************************";
	cout << endl;
	
	return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



