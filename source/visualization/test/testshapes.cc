// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testshapes.cc,v 1.1 1999-01-08 16:33:59 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  Jan 1996

#include "BuildShapes.hh"

#include "G4VisManager.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UIWo.hh"

#include "G4Initializer.hh"
#include "G4RunManager.hh"
#include "G4GeometryManager.hh"

// Solve templates.
// Code in g4templates.hh controlled by the macro G4_SOLVE_TEMPLATES
#include "g4templates.hh"

main (int argc,char** argv) {

  // Set up User Interface.
  G4UImanager* UI = G4UImanager::GetUIpointer ();

  // Choose (G)UI.
#ifdef Wo
  if(argc==2)
    {
      G4UIterminal* interactor = new G4UIterminal;
      UI -> SetSession (interactor);
    }
  else
    {
      G4UIWo*       interactor = new G4UIWo (argc,argv);
      UI -> SetSession (interactor);
    }
#else
  G4UIterminal* interactor = new G4UIterminal;
  UI -> SetSession (interactor);
#endif

  G4VPhysicalVolume* pBox  = BuildBox ();       // 0
  G4VPhysicalVolume* pCylinder = BuildCylinder();   // 1
  G4VPhysicalVolume* pTubs = BuildTubs();       // 2
  G4VPhysicalVolume* pCons = BuildCons();       // 3
  G4VPhysicalVolume* pTrd  = BuildTrd ();       // 4
  G4VPhysicalVolume* pTrap = BuildTrap();       // 5
  G4VPhysicalVolume* pSphereFull = BuildSphereFull(); // 6
  G4VPhysicalVolume* pSphereSeg  = BuildSphereSeg (); // 7
  G4VPhysicalVolume* pPara = BuildPara();       // 8
  G4VPhysicalVolume* pPCon = BuildPCon();       // 9
  G4VPhysicalVolume* pPGon = BuildPGon();       // 10
  G4VPhysicalVolume* pforcedWireframeBox = BuildForcedWireframeBox(); // 11

  G4VPhysicalVolume* pGeometries [20];
  pGeometries [0] = pBox;
  pGeometries [1] = pCylinder;
  pGeometries [2] = pTubs;
  pGeometries [3] = pCons;
  pGeometries [4] = pTrd;
  pGeometries [5] = pTrap;
  pGeometries [6] = pSphereFull ;
  pGeometries [7] = pSphereSeg ;
  pGeometries [8] = pPara ;
  pGeometries [9] = pPCon ;
  pGeometries [10] = pPGon ;
  pGeometries [11] = pforcedWireframeBox ;

  if (!pGeometries [0]) { return 0 ; }

  // Initializer, etc.
  G4Initializer * initializer = new G4Initializer;
  testshapesDetectorConstruction* pDetectorConstruction
    = new testshapesDetectorConstruction;
  pDetectorConstruction -> SetDetector (pGeometries [0]);
  initializer -> SetUserInitialization (pDetectorConstruction);
  initializer -> Initialize();

  // Instantiate Vis Manager.
  G4VisManager* pVMan = G4VisManager::GetInstance ();

  // Tie streams together so output is flushed before any input.
  cin.tie (&G4cout);
  // Main loop.
  G4int iGeom;
  do {
    // Choose a shape.
    do {
      G4cout << "Choose a shape ( < 0 to quit ):\n"
	   << "0)  Box\n"
	   << "1)  Cylinder\n"
	   << "2)  Tubs\n"
	   << "3)  Cons\n"
	   << "4)  Trd\n"
	   << "5)  Trap\n"
	   << "6)  Sphere (full)\n"
	   << "7)  Sphere (segment)\n"
	   << "8)  Pararellepiped\n"
	   << "9)  PCon (Not implemented yet)\n"
	   << "10) PGon (Not implemented yet)\n"
	   << "11) Forced wireframe box\n"
	   << endl;
      G4cout << "Enter choice (< 0 to exit): ";
      cin >> iGeom;  // (Protect against invalid input!!!!!!!!!!!!!!!!!!!)
      while (cin.get () != '\n');  // Get newline character after cin >> iGS;
      
    } while (iGeom > 11);
    
    if (iGeom >= 0) {

      if (pGeometries [iGeom]) {

	pDetectorConstruction -> SetDetector (pGeometries [iGeom]);
	initializer -> Initialize();

	pVMan -> AddToCurrentSceneData (pGeometries [iGeom]);

	// These lines temporarily removed from G4RunManager, so they're here..
	G4cout << "Start closing geometry." << endl;
	G4GeometryManager::GetInstance()->CloseGeometry();
	G4cout << "Geometry closed. Start event processing." << endl;

	// Enter UI command interpreter.
	UI -> Interact ("G4Vis> ");
      }
      else {
	G4cout << "Null geometry pointer."
	  "\nMaybe you selected a not-implemented shape!"
	     << endl;
      }
    }
  } while ( iGeom >= 0 );
  
  // Ensure VisManager destructor is called.
  delete G4VisManager::GetInstance ();
  
  return 0;
}
