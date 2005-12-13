//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// Created by Joanna Weng 26.11.2004
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
#include "ExGflashDetectorConstruction.hh"
#include "ExGflashPhysicsList.hh"
#include "ExGflashPrimaryGeneratorAction.hh"
#include "ExGflashEventAction.hh"
#include "ExGflashRunAction.hh"

#ifdef G4VIS_USE
#include "ExGflashVisManager.hh"
#endif

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
	runManager->SetUserInitialization(new ExGflashDetectorConstruction);
	cout<<"# GFlash Example: Physics list"<<endl;
	runManager->SetUserInitialization(new ExGflashPhysicsList);
	cout<<"# GFlash Example: Primary Generator"<<endl;
	runManager->SetUserAction(new ExGflashPrimaryGeneratorAction);
	cout<<"# GFlash Example: User Action Classes"<<endl;
	runManager->SetUserAction(new ExGflashEventAction);
	runManager->SetUserAction(new ExGflashRunAction);
	
	#ifdef G4VIS_USE
	G4VisManager* visManager = new ExGflashVisManager;
	visManager->Initialize();
	#endif
	
	G4UImanager* UI = G4UImanager::GetUIpointer();
	UI->ApplyCommand("/run/verbose 0");
	runManager->Initialize();
	UI->ApplyCommand("/Step/Verbose 0");
	
	if (argc==1)   // Define UI terminal for interactive mode  
	{ 
		G4UIsession * session = new G4UIterminal(new G4UItcsh);	
		UI->ApplyCommand("/control/execute vis.mac");    
		session->SessionStart();
		//delete session;
	}
	else           // Batch mode
	{ 
		G4String s=*(argv+1);
		UI->ApplyCommand("/control/execute "+s);
	}
	
	#ifdef G4VIS_USE
	delete visManager;
	#endif	
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



