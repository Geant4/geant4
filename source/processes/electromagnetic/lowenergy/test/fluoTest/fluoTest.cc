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
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif
#include "Randomize.hh"
#ifdef G4VIS_USE
#include "fluoTestVisManager.hh"
#endif
#include "fluoTestDetectorConstruction.hh"
#include "fluoTestPhysicsList.hh"
#include "fluoTestPrimaryGeneratorAction.hh"
#include "fluoTestRunAction.hh"
#include "fluoTestEventAction.hh"
#include "fluoTestSteppingAction.hh"
#include "fluoTestSteppingVerbose.hh"

#ifdef G4ANALYSIS_USE
#include "fluoTestAnalysisManager.hh"
#endif

/* This global file is used to store relevant data for
   analysis with external tools */
G4std::ofstream outFile;


int main(int argc,char** argv) {

  // choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);

  //fluoTest Verbose output class
  G4VSteppingVerbose::SetInstance(new fluoTestSteppingVerbose);
     
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  fluoTestDetectorConstruction* detector = new fluoTestDetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new fluoTestPhysicsList);
  
 G4UIsession* session=0;
  
  if (argc==1)   // Define UI session for interactive mode.
{
      // G4UIterminal is a (dumb) terminal.
#ifdef G4UI_USE_XM
      session = new G4UIXm(argc,argv);
#else           
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);
#else
      session = new G4UIterminal();
#endif
#endif
    }
  
#ifdef G4VIS_USE
  //visualization manager
   G4VisManager* visManager = new fluoTestVisManager;
  visManager->Initialize();
#endif
  
#ifdef G4ANALYSIS_USE
  // Creation of the analysis manager
  fluoTestAnalysisManager* analysisMgr = new fluoTestAnalysisManager(detector);
#endif
  
 // Set optional user action classes
#ifdef G4ANALYSIS_USE
  fluoTestEventAction* eventAction = 
    new fluoTestEventAction(analysisMgr);
  fluoTestRunAction* runAction =
    new fluoTestRunAction(analysisMgr);
  fluoTestSteppingAction* stepAction = 
    new fluoTestSteppingAction(detector,analysisMgr); 
 //#else 
  //fluoTestEventAction* eventAction = new fluoTestEventAction();
  //fluoTestRunAction* runAction = new fluoTestRunAction();
#endif 

// set user action classes
  runManager->SetUserAction(new fluoTestPrimaryGeneratorAction(detector));
  
  runManager->SetUserAction(eventAction); 
  runManager->SetUserAction(runAction);
#ifdef G4ANALYSIS_USE 
runManager->SetUserAction(stepAction);
#else
runManager->SetUserAction(new fluoTestSteppingAction);  
#endif
 //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();
 
  if (session)   // Define UI session for interactive mode.
    {
      // G4UIterminal is a (dumb) terminal.
      UI->ApplyCommand("/control/execute initInter.mac");    
#ifdef G4UI_USE_XM
      // Customize the G4UIXm menubar with a macro file :
      UI->ApplyCommand("/control/execute gui.mac");
#endif
      session->SessionStart();
      delete session;
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
 UI->ApplyCommand(command+fileName);
    }

  // job termination
#ifdef G4VIS_USE
   delete visManager;
  #endif
#ifdef G4ANALYSIS_USE
  delete analysisMgr;  
 #endif
delete runManager;
  return 0;
}

 

