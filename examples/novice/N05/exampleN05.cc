// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: exampleN05.cc,v 1.2 1999-05-28 12:55:34 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleN05
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
// --------------------------------------------------------------
// Comments
//
// * Example of a main program making use of parameterisation
//   ie "Fast Simulation".
//-------------------------------------------------------------------

//------------------------------
// Run, Generator etc.. actions:
//------------------------------
#include "ExN05RunAction.hh"
#include "ExN05PrimaryGeneratorAction.hh"
#include "ExN05EventAction.hh"
#include "ExN05SteppingAction.hh"

//---------------------------
// Parameterisation manager:
//---------------------------
#include "G4GlobalFastSimulationManager.hh"

//--------------------
// Detector:
//--------------------
#include "ExN05DetectorConstruction.hh"

//----------------------------------
// ExN05PhysicsList makes use of the
// G4ParameterisationManagerProcess
//----------------------------------
#include "ExN05PhysicsList.hh"

#ifdef GNU_GCC
  #include "ExN05CalorimeterHit.hh"
  template class RWTValVector<ExN05CalorimeterHit>;
  template class G4Allocator<ExN05CalorimeterHit>;
#endif

#include "G4UIterminal.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh"

#ifdef G4VIS_USE
#include "MyVisManager.cc"
#endif

#include "G4ios.hh"


int main(int argc, char** argv) {

  //-------------------------------
  // Initialization of Run manager
  //-------------------------------
  G4cout << "RunManager construction starting...." << endl;
  G4RunManager * runManager = new G4RunManager;

  // Detector geometry
  ExN05DetectorConstruction* Detector = new ExN05DetectorConstruction();
  runManager->SetUserInitialization(Detector);

  // PhysicsList: (including G4FastSimulationManagerProcess)
  runManager->SetUserInitialization(new ExN05PhysicsList);

  // UserAction classes.
  runManager->SetUserAction(new ExN05RunAction);
  runManager->SetUserAction(new ExN05PrimaryGeneratorAction);
  runManager->SetUserAction(new ExN05EventAction);
  runManager->SetUserAction(new ExN05SteppingAction);

  // Inizialize Run manager
  runManager->Initialize();

  // Close the "fast simulation": will
  // trigger the ghost geomtries construction:
  G4GlobalFastSimulationManager::GetGlobalFastSimulationManager()->
    CloseFastSimulation();

    
  //----------------
  // Visualization:
  //----------------
#ifdef G4VIS_USE
  G4cout << "Instantiating MyVisManager......." << endl;
  G4VisManager* visManager = new MyVisManager;
  visManager -> Initialize ();
#endif

  //Setup commandes:
  G4UImanager * UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/Step/Verbose 0");
  UI->ApplyCommand("/tracking/Verbose 1");
  UI->ApplyCommand("/gun/particle e-");
  UI->ApplyCommand("/gun/energy 100 MeV");
  UI->ApplyCommand("/gun/direction 0 0 1");
  UI->ApplyCommand("/gun/position 0 0 0");
  UI->ApplyCommand("/gun/direction 0 .3 1.");

  if(argc==1)
  {
    //--------------------------
    // Define (G)UI
    //--------------------------
    // G4UIterminal is a (dumb) terminal.
    G4UIsession* session = new G4UIterminal;
    session->SessionStart();
    delete session;
  }
  else
  {
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

  delete runManager;
#ifdef G4VIS_USE
  delete visManager;
#endif
  return EXIT_SUCCESS;
}




