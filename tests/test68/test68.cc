#include "G4RunManager.hh" 
#include "G4UImanager.hh" 
#include "Tst68DetectorConstruction.hh" 
#include "FTFP_BERT.hh"
#include "Tst68PrimaryGeneratorAction.hh" 
#include "Tst68EventAction.hh" 
#include "Tst68RunAction.hh" 
#include "Tst68TrackingAction.hh" 
#include "Tst68StackingAction.hh" 
#include "CLHEP/Random/Ranlux64Engine.h" 
#include "CLHEP/Random/MTwistEngine.h" 
#include <ctime>


int main(int argc,char** argv) { 

  CLHEP::Ranlux64Engine defaultEngine( 1234567, 4 ); 
 
  CLHEP::HepRandom::setTheEngine( &defaultEngine ); 
  G4int seed = std::time( NULL ); 
  CLHEP::HepRandom::setTheSeed( seed ); 
  G4cout << G4endl 
         << " ===================================================== " << G4endl 
         << " Initial seed = " << seed << G4endl 
	 << " ===================================================== " << G4endl 
	 << G4endl; 
  G4RunManager* runManager = new G4RunManager; 
  runManager->SetUserInitialization( new Tst68DetectorConstruction ); 

  FTFP_BERT *thePL = new FTFP_BERT;

  //thePL->SetDefaultCutValue( 0.020 *mm ); // 20 microns 
  runManager->SetUserInitialization( thePL ); 
  runManager->SetUserAction( new Tst68PrimaryGeneratorAction ); 
  runManager->SetUserAction( new Tst68RunAction ); 
  runManager->SetUserAction( new Tst68EventAction ); 
  runManager->SetUserAction( new Tst68TrackingAction ); 
  runManager->SetUserAction( new Tst68StackingAction ); 

  runManager->Initialize(); 
  G4UImanager* UI = G4UImanager::GetUIpointer(); 
  if ( argc==1 ) {   // Define UI session for interactive mode. 
  } else {   // Batch mode 
    G4String command = "/control/execute "; 
    G4String fileName = argv[1]; 
    UI->ApplyCommand(command+fileName); 
  } 

  // job termination 
  delete runManager; 
  return 0; 
} 
