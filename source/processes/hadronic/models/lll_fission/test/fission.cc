//******************************************************************************
// fission.cc  GEANT4 user application for testing handling of
//             low-energy neutron induced fission
//
// fission                 # start interactive version
// fission test.mac        # read test.mac for geant commands
// fission -i test.mac     # read test.mac for geant commans and continue interactively
//
// 1.00 JMV, LLNL, JAN-2007:  First version.
//******************************************************************************
//

// misc includes
//
#include <fstream>
#include <math.h>
#include "G4ios.hh"

// package includes
//
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"

// geant4 includes
//
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

// visualization
//
#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

//------------------------------------------------------------------------------
int main(int argc,char** argv) {

  G4String filename;
  bool interactive=false;
  G4String interactiveFlag = "-i";

 //....Last option is filename (unless it is -i)
  if(argc>1 && G4String(argv[argc-1]) != interactiveFlag) 
    filename = G4String(argv[argc-1]); 

 //....Make session interactive if "-i"
  if(argc==1 || G4String(argv[1]) == interactiveFlag) interactive = true;


  // Run manager
  //------------
  G4RunManager * theRunManager = new G4RunManager;

  // UserInitialization classes
  //---------------------------
  DetectorConstruction* theDetector = new DetectorConstruction;
  PhysicsList* thePhysicsList = new PhysicsList;

  theRunManager->SetUserInitialization(theDetector);
  theRunManager->SetUserInitialization(thePhysicsList);

  // UserAction classes
  //-------------------
  theRunManager->SetUserAction(new PrimaryGeneratorAction());
  
#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  //-----------------------------------------
  G4VisManager* theVisManager = new G4VisExecutive;
  theVisManager->Initialize();
#endif
   
  // Initialize G4 kernel
  //---------------------
  theRunManager->Initialize();

  // User interactions
  //------------------
  if(filename) 
    G4UImanager::GetUIpointer()->ApplyCommand("/control/execute "+filename);

  //....Start an interactive session
  if(interactive){
    G4cout << "\nType 'exit' to end the program.\n";
    G4UIterminal(new G4UItcsh).SessionStart();
  }

#ifdef G4VIS_USE
  delete theVisManager;
#endif
  delete theRunManager;
  return EXIT_SUCCESS;
}
