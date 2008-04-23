#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIsession.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "CLHEP/Random/RanluxEngine.h" 

#include "MyDetectorConstruction.hh" 
#include "MyPrimaryGeneratorAction.hh" 
#include "MyEventAction.hh" 
#include "LHEP.hh" 
#include "QGSP.hh" 
#include "QGSP_BERT.hh" 
#include "QGSP_BIC.hh" 
#include "QGSC.hh" 
#include "FTFP.hh" 
#include "FTFC.hh" 
#include "QGSP_BERT_HP.hh" 
#include "QGSP_EMV.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "G4GeometryManager.hh"


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

  // Added to test Gabriele's tolerance.
  // The argument is the world's maximum length, which is printed out
  // after building the GDML geometry (see MyDetectorConstruction.cc).
  // G4GeometryManager::GetInstance()->SetWorldMaximumExtent( 54.0*m );

#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  runManager->SetUserInitialization( new MyDetectorConstruction ); 

  LHEP *thePL = new LHEP; //***LOOKHERE*** SELECT PHYSICS LIST.

  //thePL->SetDefaultCutValue( 0.020 *mm ); // 20 microns 
  runManager->SetUserInitialization( thePL ); 
  runManager->SetUserAction( new MyPrimaryGeneratorAction ); 
  runManager->SetUserAction( new MyEventAction ); 

  runManager->Initialize(); 

  // get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();
                                                                                
  if ( argc==1 ) {   // Define UI session for interactive mode.
    G4UIsession * session = new G4UIterminal(new G4UItcsh);
    UI->ApplyCommand("/control/execute vis.mac");
    G4cout << "Now, please, apply beamOn command..." << G4endl;
    session->SessionStart();
    delete session;
  } else {   // Batch mode 
    G4String command = "/control/execute "; 
    G4String fileName = argv[1]; 
    UI->ApplyCommand(command+fileName); 
  } 

  G4cout << G4endl 
	 << " ===================================================== " << G4endl 
         << " Final random number = " 
         << CLHEP::HepRandom::getTheEngine()->flat() << G4endl 
	 << " ===================================================== " << G4endl 
         << G4endl; 

  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0; 
} 
