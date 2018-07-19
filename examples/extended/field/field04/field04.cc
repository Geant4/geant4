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
// $Id: field04.cc 110170 2018-05-17 07:28:42Z gcosmo $
//
/// \file field/field04/field04.cc
/// \brief Main program of the field/field04 example
//
//
// --------------------------------------------------------------
//
//      GEANT 4 - Example F04
//
// --------------------------------------------------------------
// Comments
//
//
// --------------------------------------------------------------

#ifndef WIN32
#include <unistd.h>
#endif

#include "G4Types.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "F04SteppingVerbose.hh"
#include "G4RunManager.hh"
#endif

#include "F04PhysicsList.hh"
#include "F04DetectorConstruction.hh"

#include "F04ActionInitialization.hh"

#include "G4UImanager.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

// argc holds the number of arguments (including the name) on the command line
// -> it is ONE when only the name is  given !!!
// argv[0] is always the name of the program
// argv[1] points to the first argument, and so on

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Instantiate G4UIExecutive if there are no arguments (interactive mode)
  G4UIExecutive* ui = nullptr;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  G4int myseed = 1234;
  if (argc  > 2) myseed = atoi(argv[argc-1]);

  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
#else
  G4VSteppingVerbose::SetInstance(new F04SteppingVerbose);
  G4RunManager * runManager = new G4RunManager;
#endif

  G4Random::setTheSeed(myseed);

  G4String physicsList = "QGSP_BERT";

#ifndef WIN32
  G4int c = 0;
  while ((c=getopt(argc,argv,"p")) != -1)
  {
     switch (c)
     {
       case 'p':
         physicsList = optarg;
         G4cout << "Physics List used is " <<  physicsList << G4endl;
         break;
       case ':':       /* -p without operand */
         fprintf(stderr,"Option -%c requires an operand\n", optopt);
         break;
       case '?':
         fprintf(stderr,"Unrecognised option: -%c\n", optopt);
     }
  }
#endif

  // Set mandatory initialization classes
  //
  // Detector construction
  F04DetectorConstruction* detector = new F04DetectorConstruction();
  runManager->SetUserInitialization(detector);
  // Physics list
  runManager->SetUserInitialization(new F04PhysicsList(physicsList));
  // User action initialization
  runManager->SetUserInitialization(new F04ActionInitialization(detector));

  // Initialize G4 kernel
  //
  //runManager->Initialize();

  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

#ifndef WIN32
  G4int optmax = argc;
  if (argc > 2)  { optmax = optmax-1; }

  if (optind < optmax)
  {
     G4String command = "/control/execute ";
     for ( ; optind < optmax; optind++)
     {
         G4String macroFilename = argv[optind];
         UImanager->ApplyCommand(command+macroFilename);
     }
  }
#else  // Simple UI for Windows runs, no possibility of additional arguments
  if (!ui)   // batch mode
  {
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UImanager->ApplyCommand(command+fileName);
  }
#endif
  else
  {
//     UImanager->ApplyCommand("/control/execute vis.mac");
     G4cout << "At the prompt, issue commands to set up detector & field, then:"
         << G4endl;
     G4cout << "/run/initialize" << G4endl;
     G4cout << "Then if you want a viewer:"<< G4endl;
     G4cout << "/control/execute vis.mac" << G4endl;
     G4cout << "Then: " << G4endl;
     G4cout << "/run/beamOn … etc." << G4endl;
     if (ui->IsGUI())
        UImanager->ApplyCommand("/control/execute gui.mac");
     ui->SessionStart();
     delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

  delete visManager;
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
