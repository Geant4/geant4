// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTel.cc,v 1.2 2000-11-15 20:27:38 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 main program
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTel example main program ------
//           by F.Longo, R.Giannitrapani & G.Santin (13 nov 2000)
//
// ************************************************************

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif


#ifdef G4VIS_USE
#include "GammaRayTelVisManager.hh"
#endif

#include "GammaRayTelDetectorConstruction.hh"
#include "GammaRayTelPhysicsList.hh"
#include "GammaRayTelPrimaryGeneratorAction.hh"
#include "GammaRayTelRunAction.hh"
#include "GammaRayTelEventAction.hh"

#ifdef G4HIS_USE_AIDA
#include "GammaRayTelHistogram.hh"
#endif

#include "g4std/vector"

// This global file is used to store relevant data for
// analysis with a separate program

G4std::ofstream outFile("tracks.dat");

G4bool drawEvent;
G4std::vector<G4String> EnteringParticles;
G4std::vector<G4double> EnteringEnergy;
G4std::vector<G4ThreeVector> EnteringDirection;

int main(int argc, char** argv)
{
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;
  
  // set mandatory initialization classes
  GammaRayTelDetectorConstruction* detector = new GammaRayTelDetectorConstruction;
  runManager->SetUserInitialization(detector);
  
  runManager->SetUserInitialization(new GammaRayTelPhysicsList);

  
#ifdef G4HIS_USE_AIDA
  // Creation of the histogram manager
  GammaRayTelHistogram histoMgr(detector);
#endif
  
  G4UIsession* session=0;
  
#ifdef G4UI_USE_XM
  session = new G4UIXm(argc,argv);
#else
  session = new G4UIterminal;
#endif
  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new GammaRayTelVisManager;
  visManager->Initialize();
#endif

  // set mandatory user action classes
  runManager->SetUserAction(new GammaRayTelPrimaryGeneratorAction(detector));
  runManager->SetUserAction(new GammaRayTelRunAction);
  
  // set optional user action classes
#ifdef G4HIS_USE_AIDA
  GammaRayTelEventAction* eventAction = new GammaRayTelEventAction(&histoMgr);
#else 
  GammaRayTelEventAction* eventAction = new GammaRayTelEventAction();
#endif
  
  runManager->SetUserAction(eventAction);
  
  // Initialize G4 kernel
  runManager->Initialize();
  

  // get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();
  
  if (session) 
    {
      // prerunGammaRayTel.mac is loaded by default
      // if macro file i exists, it is passed as the argument
      // of the executable

      if(argc>1)
	{
	  G4String command = "/control/execute ";
	  for (int i=2; i<=argc; i++) 
	    {
	      G4String macroFileName = argv[i-1];
	      UI->ApplyCommand(command+macroFileName);
	    }
	}
      else  UI->ApplyCommand("/control/execute prerunGammaRayTel.mac");
      session->SessionStart();
      delete session;
    }

  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  return 0;
}








