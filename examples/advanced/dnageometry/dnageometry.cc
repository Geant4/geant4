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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// and the DNA geometry given in the Geom_DNA example 
// shall cite the following Geant4-DNA collaboration publications:
// [1] NIM B 298 (2013) 47-54
// [2] Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4RunManager.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4VIS_USE
  #include "G4VisExecutive.hh"
#endif

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

int main(int argc,char** argv) 
{
  // Choose the Random engine
  
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // Set mandatory user initialization classes
  DetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization(detector);
  
  runManager->SetUserInitialization(new PhysicsList);

  // Set mandatory user action classes
  runManager->SetUserAction(new PrimaryGeneratorAction(detector));
  PrimaryGeneratorAction* primary = new PrimaryGeneratorAction(detector);

  HistoManager*  histo = new HistoManager();

  // Set optional user action classes
  RunAction* RunAct = new RunAction(detector,histo);
  runManager->SetUserAction(RunAct);

  runManager->SetUserAction(new SteppingAction(RunAct,detector,primary,histo));
  
#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
   
  // Initialize G4 kernel
  runManager->Initialize();
    
  remove ("geom_dna.root");  
    
  // Get the pointer to the User Interface manager 
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  

  if (argc==1)   // Define UI session for interactive mode.
  { 
#ifdef _WIN32
    G4UIsession * session = new G4UIterminal();
#else
    G4UIsession * session = new G4UIterminal(new G4UItcsh);
#endif
    UImanager->ApplyCommand("/control/execute dnageometry.mac");    
    session->SessionStart();
    delete session;
  }
  else           // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }

#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;

  return 0;
}

