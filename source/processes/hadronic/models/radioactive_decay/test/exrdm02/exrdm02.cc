//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: exrdm02.cc,v 1.1 2003-10-08 16:31:48 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - GGEVis.cc
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#ifdef G4UI_USE_GAG
#include "G4UIGAG.hh"
#endif
#include "G4UIterminal.hh"
#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif 

#include "exrdm02GeometryConstruction.hh"
#include "exrdm02PhysicsList.hh"
#include "exrdm02EventAction.hh"
#include "exrdm02RunAction.hh"
#include "exrdm02SteppingAction.hh"
#include "exrdm02PrimaryGeneratorAction.hh"
#include "Randomize.hh"

#include <vector>
//G4String filename;
G4bool drawEvent;
std::vector<G4String> Particles;
std::vector<G4double> Energies;
std::vector<G4double> Weights;
std::vector<G4double> Times;

#ifdef G4VIS_USE
#include "exrdm02VisManager.hh"
#endif

int main(int argc,char** argv)
{

  //  G4cout << " The results file name = "<<G4endl ;
  //G4cin >> filename;

  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // set mandatory initialization classes

  exrdm02GeometryConstruction* Geometry = new exrdm02GeometryConstruction;
  runManager->SetUserInitialization(Geometry);
  runManager->SetUserInitialization(new exrdm02PhysicsList);

#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new exrdm02VisManager;
  visManager->Initialize();
#endif

  // set mandatory user action class
  runManager->SetUserAction(new exrdm02PrimaryGeneratorAction);
  runManager->SetUserAction(new exrdm02RunAction);
  runManager->SetUserAction(new exrdm02EventAction);
  runManager->SetUserAction(new exrdm02SteppingAction);
  // Initialize G4 kernel
  runManager->Initialize();

  // get the pointer to the User Interface manager 
    
      G4UImanager* UI = G4UImanager::GetUIpointer();  
      UI->ApplyCommand("/run/verbose 1");
      UI->ApplyCommand("/event/verbose 2");
      UI->ApplyCommand("/tracking/verbose 1");           
      if (argc==1) {
	//     	G4UIsession * session = new G4UIGAG;
     	G4UIsession * session = new G4UIterminal;
      	session->SessionStart();
      	delete session;
  }
      else {
 
  //
  //
  // Create a pointer to the user interface manager.
  //
	//      G4UImanager *UI = G4UImanager::GetUIpointer();
      G4String command = "/control/execute ";
      for (int i=2; i<=argc; i++) {
      G4String macroFileName = argv[i-1];
      UI->ApplyCommand(command+macroFileName);
       }

  }                                  

  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;

  return 0;
}








