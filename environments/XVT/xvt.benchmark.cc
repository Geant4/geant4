// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: xvt.benchmark.cc,v 1.3 1999-12-15 14:48:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Updated for XVT GUI by Simon Prior 12/08/97

#include "MyDetectorConstruction.cc"
#include "MyRunAction.cc"
#include "MyPrimaryGeneratorAction.cc"
#include "MyEventAction.cc"
#include "MySteppingAction.cc"

#include "G4UIxvt.hh"
#include "G4Initializer.hh"
#include "G4RunManager.hh"

#ifdef G4VIS_USE
  #include "MyVisManager.cc"
#endif

#include "G4ios.hh"

// Solve templates.
#ifdef GNU_GCC
  #include "MyCalorimeterHit.hh"
  template class G4RWTValVector<MyCalorimeterHit>;
  template class G4Allocator<MyCalorimeterHit>;
#endif
#include "g4templates.hh"

main() {


  // Initializer
  G4Initializer * initializer = new G4Initializer;

  // Detector geometry
  initializer->SetUserInitialization(new MyDetectorConstruction);

  // Initialization
  initializer->Initialize();

  // Run manager
  G4cout << "RunManager construction starting...." << G4endl;
  G4RunManager * runManager = new G4RunManager;

  // UserAction classes.
  runManager->SetUserAction(new MyRunAction);
  runManager->SetUserAction(new MyPrimaryGeneratorAction);
  runManager->SetUserAction(new MyEventAction);
  runManager->SetUserAction(new MySteppingAction);

#ifdef G4VIS_USE
  G4cout << "Instantiating MyVisManager......." << G4endl;
  G4VisManager* visManager = new MyVisManager;
  visManager -> Initialize ();
#endif

  // Event loop
  // Define (G)UI XVT
  G4UIsession * session = new G4UIxvt;
  session->SessionStart();

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete session;
  delete initializer;
  delete runManager;
  return 0;
}

