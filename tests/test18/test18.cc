//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$
//
// --------------------------------------------------------------
//      GEANT 4 - test18.cc
// --------------------------------------------------------------

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "Tst18ActionInitialization.hh"
#include "Tst18GeometryConstruction.hh"
#include "Tst18PhysicsList.hh"
#include "Randomize.hh"

#include <vector>

std::vector<G4String> Particles;
std::vector<G4double> Energies;
std::vector<G4double> Weights;
std::vector<G4double> Times;

int main(int argc,char** argv)
{
  // Run manager
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(1);
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  // set mandatory initialization classes

  runManager->SetUserInitialization(new Tst18GeometryConstruction);
  runManager->SetUserInitialization(new Tst18PhysicsList);
  runManager->SetUserInitialization(new Tst18ActionInitialization);

  runManager->Initialize();

  // get the pointer to the User Interface manager 
    
  G4UImanager* UI = G4UImanager::GetUIpointer();  
  UI->ApplyCommand("/run/verbose 0");
  UI->ApplyCommand("/event/verbose 0");
  UI->ApplyCommand("/tracking/verbose 0");           
  UI->ApplyCommand("/process/verbose 0");           
  if (argc==1) {
    G4UIsession * session = new G4UIterminal;
    session->SessionStart();
    delete session;
  }
  else {
 
    // Create a pointer to the user interface manager.
    G4String command = "/control/execute ";
    for (int i=2; i<=argc; i++) {
       G4String macroFileName = argv[i-1];
       UI->ApplyCommand(command+macroFileName);
    }
  }                                  

  // job termination
  delete runManager;
  return 0;
}
