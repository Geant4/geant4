// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTel.cc                            main file *     
// * -------                                                            *
// *                                                                    *
// * Version:           0.4                                             *
// * Date:              06/11/00                                        *
// * Author:            R Nartallo                                      *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
//
// HISTORY
// -------
//
// The development of this advanced example is based on earlier work 
// carried out by a team of Geant4 collaborators  to simulate the Chandra
// and XMM X-ray observatories. The authors involved in those models are
// J Apostolakis, P Arce, S Giani, F Lei, R Nartallo, S Magni,
// P Truscott, L Urban
//
// **********************************************************************
// 
// CHANGE HISTORY
// --------------
//
// 06.11.2000 R.Nartallo
// - First implementation of X-ray Telescope advanced example.
// - Based on Chandra and XMM models
// - Lines for using GAG and the histogram manager are commented out.
//
//
// **********************************************************************

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "XrayTelDetectorConstruction.hh"
#include "XrayTelPhysicsList.hh"
#include "XrayTelEventAction.hh"
#include "XrayTelRunAction.hh"
#include "XrayTelSteppingAction.hh"
#include "XrayTelPrimaryGeneratorAction.hh"

//#include "G4UIterminal.hh"
#ifdef G4UI_USE_GAG
#include "G4UIGAG.hh"
#endif
#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif
#ifdef G4VIS_USE
#include "XrayTelVisManager.hh"
#endif                                                                                    
#ifdef G4ANALYSIS_USE
#include "XrayTelAnalysisManager.hh"
#endif

int main( int argc, char** argv )
{
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  XrayTelDetectorConstruction* telescope = new XrayTelDetectorConstruction;
  runManager->SetUserInitialization( telescope );
  runManager->SetUserInitialization(new XrayTelPhysicsList);

  // set mandatory user action class
  runManager->SetUserAction(new XrayTelPrimaryGeneratorAction( telescope ));

  // create a manager for the analysis
  // The analysis system is given through an environment variable :
#ifdef G4ANALYSIS_USE
  char* s = getenv("G4ANALYSIS_SYSTEM");
  XrayTelAnalysisManager* analysisManager = new XrayTelAnalysisManager(s?s:"");
  runManager->SetUserAction(new XrayTelRunAction(analysisManager));
  runManager->SetUserAction(new XrayTelEventAction(analysisManager));
  runManager->SetUserAction(new XrayTelSteppingAction(analysisManager));
#else
  runManager->SetUserAction(new XrayTelRunAction());
  runManager->SetUserAction(new XrayTelEventAction());
  runManager->SetUserAction(new XrayTelSteppingAction());
#endif

   
  // visualization manager
#ifdef G4VIS_USE
  G4VisManager* visManager = new XrayTelVisManager;
  visManager->Initialize();    
#endif

  //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4UImanager *UI = G4UImanager::GetUIpointer();  
  if ( argc==1 ){
    G4UIsession* session = 0;
#ifdef G4UI_USE_XM
    session = new G4UIXm(argc,argv);
    // Customize the G4UIXm menubar with a macro file :
    UI->ApplyCommand("/control/execute gui.mac");
#else
    //session = new G4UIterminal;
    session = new G4UIGAG;
#endif
    session->SessionStart();
    delete session;
  } else {
    // Create a pointer to the User Interface manager 
    G4String command = "/control/execute ";
    for (int i=2; i<=argc; i++) {
      G4String macroFileName = argv[i-1];
      UI->ApplyCommand(command+macroFileName);
    }
  }                                  

  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
#ifdef G4ANALYSIS_USE
  delete analysisManager;
#endif
  delete runManager;
  return 0;
}

