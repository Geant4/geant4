#include "G4RunManager.hh" 
#include "G4UImanager.hh" 
#include "CLHEP/Random/RanluxEngine.h" 
#include "StatAccepTestDetectorConstruction.hh" 
#include "StatAccepTestPrimaryGeneratorAction.hh" 
#include "StatAccepTestEventAction.hh" 
#include "StatAccepTestRunAction.hh" 
#include "StatAccepTestTrackingAction.hh" 
#include "StatAccepTestStackingAction.hh" 
#include "StatAccepTestAnalysis.hh" 
#include "G4UIterminal.hh" 
#include <cstdlib>
#include "time.h"

//***LOOKHERE***
// Uncomment these includes files only if you use a Geant4 version before 9.2 .
// #include "LHEP.hh"
// #include "QGSP.hh"
// #include "QGSP_BERT.hh"
// #include "QGSP_BIC.hh"
// #include "QGSC.hh"
// #include "QGSP_EMV.hh"
// #include "FTFP.hh"
// #include "FTFC.hh"
// Comment these includes files only if you use a Geant4 version before 9.2 .
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"

#ifdef G4UI_USE_TCSH 
  #include "G4UItcsh.hh" 
#endif 

#ifdef G4VIS_USE 
  #include "StatAccepTestVisManager.hh" 
#endif 

#ifdef G4UI_USE_XM 
#include "G4UIXm.hh" 
#endif 


int main(int argc,char** argv) { 

  CLHEP::RanluxEngine defaultEngine( 1234567, 4 ); 
  CLHEP::HepRandom::setTheEngine( &defaultEngine ); 
  G4int seed = time( NULL ); 
  CLHEP::HepRandom::setTheSeed( seed ); 
  G4cout << G4endl 
         << " ===================================================== " << G4endl 
         << " Initial seed = " << seed << G4endl 
	 << " ===================================================== " << G4endl 
	 << G4endl; 

  G4RunManager* runManager = new G4RunManager; 
  runManager->SetUserInitialization( new StatAccepTestDetectorConstruction ); 

  // ***LOOKHERE***
  // If you use a Geant4 version before 9.2, then you have to write
  // the Physics List explicitly below:
  // QGSP  *thePL = new QGSP;
  // and eventually set the production range cut with:
  // //thePL->SetDefaultCutValue( 1.0*cm );
  // You need finally the following line:
  // runManager->SetUserInitialization( thePL );
  // 
  // Starting instead with Geant4 version 9.2, you can create a
  // single executable, and then select the Physics List at run
  // time, via the PHYSLIST environmental variable.
  // For instance:
  //                             export PHYSLIST=FTFP
  // sets the FTFP Physics List.
  // If $PHYSLIST is undefined, then by default the Physics List 
  //  QGSP_BERT will be automatically selected.   
  G4PhysListFactory factory;
  G4VModularPhysicsList* physList = factory.ReferencePhysList();
  runManager->SetUserInitialization( physList );

  runManager->SetUserAction( new StatAccepTestPrimaryGeneratorAction ); 
  runManager->SetUserAction( new StatAccepTestRunAction ); 
  runManager->SetUserAction( new StatAccepTestEventAction ); 
  runManager->SetUserAction( new StatAccepTestTrackingAction ); 

  char* nameEnvironmentalForBiasing = getenv( "IS_BIASING_ACTIVE" );
  if ( nameEnvironmentalForBiasing ) {
  G4cout << " *** nameEnvironmentalForBiasing=" << nameEnvironmentalForBiasing 
	 << " ***" << G4endl << G4endl;
    runManager->SetUserAction( new StatAccepTestStackingAction ); 
  }

#ifdef flagHISTOGRAMS_ON
  StatAccepTestAnalysis::getInstance()->setIsHistogramOn( true ); 
#else
  StatAccepTestAnalysis::getInstance()->setIsHistogramOn( false );
#endif

#ifdef G4VIS_USE 
  StatAccepTestVisManager *visManager = new StatAccepTestVisManager; 
  visManager->Initialize(); 
#endif 

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
         << " Final random number = "          << CLHEP::HepRandom::getTheEngine()->flat() << G4endl 
	 << " ===================================================== " << G4endl 
         << G4endl; 
  // job termination 
#ifdef G4VIS_USE 
  delete visManager; 
#endif 
  delete runManager; 
  return 0; 
} 
