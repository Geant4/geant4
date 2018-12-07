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
//
// --------------------------------------------------------------
//      GEANT 4 - exampleN03
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------

#include "G4Types.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"

#include "G4UIExecutive.hh"

#include "FCALTestbeamSetup.hh"
#include "G4PhysListFactory.hh"
#include "FCALActionInitialization.hh"


int main(int argc,char** argv) {

  // choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
#else
    G4RunManager* runManager = new G4RunManager;
#endif

  // set mandatory initialization classes
  FCALTestbeamSetup* detector = new FCALTestbeamSetup;
  runManager->SetUserInitialization(detector);

  G4PhysListFactory factory;
  runManager->SetUserInitialization(factory.ReferencePhysList());

  runManager->SetUserInitialization(new FCALActionInitialization);

  // get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  if (argc==1)   // Define UI session for interactive mode.
     {
       G4UIExecutive* ui = new G4UIExecutive(argc, argv);
       UImanager->ApplyCommand("/control/execute prerunlArcal.mac");
       if (ui->IsGUI())
	 UImanager->ApplyCommand("/control/execute gui.mac");
       ui->SessionStart();
       delete ui;
    }
  else           // Batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }

  // job termination
  delete visManager;
  delete runManager;

  return 0;
}

