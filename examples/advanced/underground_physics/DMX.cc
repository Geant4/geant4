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
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//            Underground Advanced example main program
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// main program
// --------------------------------------------------------------

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "DMXAnalysisManager.hh"
#include "DMXDetectorConstruction.hh"
#include "DMXPhysicsList.hh"
#include "DMXActionInitializer.hh"

int main(int argc,char** argv) {

  // choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
  // Construct the default run manager
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  //runManager->SetNumberOfThreads(2); 
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  // set mandatory initialization classes
  runManager->SetUserInitialization(new DMXDetectorConstruction);
  runManager->SetUserInitialization(new DMXPhysicsList);
  runManager->SetUserInitialization(new DMXActionInitializer());

#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif


#ifdef DMXENV_GPS_USE
  G4cout << " Using GPS and not DMX gun " << G4endl;
#else
  G4cout << " Using the DMX gun " << G4endl;
#endif
    
 
  //Initialize G4 kernel
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4UImanager* UImanager = G4UImanager::GetUIpointer();  

  // Define UI session for interactive mode.
  if(argc == 1)
    {
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute initInter.mac");     
#endif
      ui->SessionStart();
      delete ui;
#endif
    }
  // Batch mode
  else
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }

  //Close-out analysis:
  // Save histograms
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();
  // Complete clean-up
  delete G4AnalysisManager::Instance();

#ifdef G4VIS_USE
  if(visManager) delete visManager;
#endif
  delete runManager;

  return 0;
}

