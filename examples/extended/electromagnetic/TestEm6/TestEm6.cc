/* hbu

cd $G4INSTALL/examples/extended/electromagnetic/TestEm6/
gmake
rm -f /tmp/hbu/*out* ; $G4WORKDIR/bin/Linux-egcs/TestEm6 run01.mac > /tmp/hbu/TestEm6.out

or

bsub -J TestEm601 -q si_1nh subm_TestEm6

or

gmake visclean ; rm -f core ; $G4WORKDIR/bin/Linux-egcs/TestEm6
/control/execute vis.mac

*/

//  TestEm6 by H.Burkhardt, April 2002

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
  #include "Em6VisManager.hh"
#endif

#include "Em6DetectorConstruction.hh"
#include "Em6PhysicsList.hh"
#include "Em6PrimaryGeneratorAction.hh"
#include "Em6RunAction.hh"
#include "Em6EventAction.hh"
#include "Em6TrackingAction.hh"
#include "Em6SteppingAction.hh"
#include "Em6SteppingVerbose.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {
 
  //choose the Random engine
  HepRandom::setTheEngine(new RanecuEngine);
  
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new Em6SteppingVerbose);
     
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  Em6DetectorConstruction* detector = new Em6DetectorConstruction;
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new Em6PhysicsList);
  
  Em6PrimaryGeneratorAction* primary = new Em6PrimaryGeneratorAction(detector);
  runManager->SetUserAction(primary);
    
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new Em6VisManager;
  visManager->Initialize();
#endif
    
  // set user action classes
  Em6RunAction* RunAct = new Em6RunAction(detector,primary);
  runManager->SetUserAction(RunAct);
  runManager->SetUserAction(new Em6EventAction   (RunAct));
  runManager->SetUserAction(new Em6TrackingAction(RunAct));
  runManager->SetUserAction(new Em6SteppingAction(detector,RunAct)); 
  
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (argc==1)   // Define UI terminal for interactive mode.
    {
      G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
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
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..... 
