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
#include "XrayFluoVisManager.hh"
#endif
#include "XrayFluoDetectorConstruction.hh"
#include "XrayFluoPhysicsList.hh"
#include "XrayFluoPrimaryGeneratorAction.hh"
#include "XrayFluoRunAction.hh"
#include "XrayFluoEventAction.hh"
#include "XrayFluoSteppingAction.hh"
#include "XrayFluoSteppingVerbose.hh"

#ifdef G4ANALYSIS_USE
#include "XrayFluoAnalysisManager.hh"
#endif

/* This global file is used to store relevant data for
   analysis with external tools */
G4std::ofstream outFile;
G4std::ofstream outFileGamma;

int main(int argc,char** argv) {

  // choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);

  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new XrayFluoSteppingVerbose);
     
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  XrayFluoDetectorConstruction* detector = new XrayFluoDetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new XrayFluoPhysicsList);
  
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
   G4VisManager* visManager = new XrayFluoVisManager;
  visManager->Initialize();
  #endif
  
#ifdef G4ANALYSIS_USE
  // Creation of the analysis manager
  XrayFluoAnalysisManager* analysisMgr = new XrayFluoAnalysisManager(detector);
#endif
  
 // Set optional user action classes
#ifdef G4ANALYSIS_USE
  XrayFluoEventAction* eventAction = 
    new XrayFluoEventAction(analysisMgr);
  XrayFluoRunAction* runAction =
    new XrayFluoRunAction(analysisMgr);
  XrayFluoSteppingAction* stepAction = 
    new XrayFluoSteppingAction(detector,analysisMgr); 

#endif 

// set user action classes

#ifdef G4ANALYSIS_USE
  runManager->SetUserAction(new XrayFluoPrimaryGeneratorAction(detector));
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(stepAction);
#else
  runManager->SetUserAction(new XrayFluoPrimaryGeneratorAction(detector));
  runManager->SetUserAction(new XrayFluoRunAction);
  runManager->SetUserAction(new XrayFluoEventAction);
  runManager->SetUserAction(new XrayFluoSteppingAction);  
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

 

