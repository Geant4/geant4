#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "StatAccepTestDetectorConstruction.hh"
#include "LHEP.hh"
#include "QGSP.hh"
#include "QGSC.hh"
#include "QGSP_BERT.hh"
#include "QGSP_BIC.hh"
#include "QGSP_GN.hh"
#include "StatAccepTestPrimaryGeneratorAction.hh"

#include "StatAccepTestEventAction.hh"
#include "StatAccepTestRunAction.hh"

#include "G4UIterminal.hh"
#ifdef G4UI_USE_TCSH
  #include "G4UItcsh.hh"
#endif
#ifdef G4VIS_USE
  #include "StatAccepTestVisManager.hh"
#endif
#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#include "CLHEP/Random/RanluxEngine.h"


int main(int argc,char** argv) {

  RanluxEngine defaultEngine( 1234567, 4 ); 
  HepRandom::setTheEngine( &defaultEngine );

  G4int seed = time( NULL );
  HepRandom::setTheSeed( seed );

  G4cout << G4endl
         << " ===================================================== " << G4endl
         << " Initial seed = " << seed << G4endl
	 << " ===================================================== " << G4endl 
	 << G4endl;

  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

#ifdef G4VIS_USE
  StatAccepTestVisManager *visManager = new StatAccepTestVisManager;
  visManager->Initialize();
#endif        

  // Set mandatory initialization classes.
  runManager->SetUserInitialization( new StatAccepTestDetectorConstruction );

  //***LOOKHERE***
  LHEP       *thePL = new LHEP;
  // QGSP       *thePL = new QGSP;
  // QGSC       *thePL = new QGSC;
  // QGSP_BERT  *thePL = new QGSP_BERT;
  // QGSP_BIC   *thePL = new QGSP_BIC;
  // QGSP_GN    *thePL = new QGSP_GN;
  // thePL->SetDefaultCutValue( 100.0*cm );
  //***endLOOKHERE***

  runManager->SetUserInitialization( thePL );

  // Set mandatory user action class.
  runManager->SetUserAction( new StatAccepTestPrimaryGeneratorAction );

  // Set optional user action classes.
  runManager->SetUserAction( new StatAccepTestRunAction );  
  runManager->SetUserAction( new StatAccepTestEventAction );
  
  // Initialize G4 kernel
  runManager->Initialize();

  G4UImanager* UI = G4UImanager::GetUIpointer();   
  if ( argc==1 ) {   // Define UI session for interactive mode.

    G4UIsession* session = 0;
#ifdef G4UI_USE_XM
    session = new G4UIXm(argc,argv);
#else
#ifdef G4UI_USE_TCSH
    session = new G4UIterminal(new G4UItcsh);      
#else
    session = new G4UIterminal();
#endif
#endif
    
#ifdef G4VIS_USE
    // Create empty scene
    G4String visCommand = "/vis/scene/create";
    UI->ApplyCommand(visCommand);
    // Choose one default viewer (you can always change it later on) 
#ifdef WIN32
    visCommand = "/vis/open VRML2FILE";
#else
    // visCommand = "/vis/open VRML2";
    visCommand = "/vis/open OGLIX";
#endif
    UI->ApplyCommand(visCommand);
    visCommand = "/vis/viewer/flush";
    UI->ApplyCommand(visCommand);
    visCommand = "/tracking/storeTrajectory 1";
    UI->ApplyCommand(visCommand);
#endif

#ifdef G4UI_USE_XM
    // Customize the G4UIXm menubar with a macro file :
    UI->ApplyCommand("/control/execute gui.g4");
#else
    G4cout << "Now, please, apply beamOn command..." << G4endl;
#endif
    
    session->SessionStart();
    delete session;
    
  } else {   // Batch mode

    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);

  }

  G4cout << G4endl 
	 << " ===================================================== " << G4endl
         << " Final random number = " << HepRandom::getTheEngine()->flat() << G4endl
	 << " ===================================================== " << G4endl 
         << G4endl;

  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  return 0;

}


