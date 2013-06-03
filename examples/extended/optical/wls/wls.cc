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
/// \file optical/wls/wls.cc
/// \brief Main program of the optical/wls example
//
//
//
//
// --------------------------------------------------------------
//      GEANT 4 - Example wls
//
// --------------------------------------------------------------
// Comments
//
//
// --------------------------------------------------------------

#ifndef WIN32
#include <unistd.h>
#endif

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#include "WLSWorkerInitialization.hh"
#else
#include "G4RunManager.hh"
#include "WLSSteppingVerbose.hh"
#endif

#include "G4UImanager.hh"

#include "Randomize.hh"

#include "WLSPhysicsList.hh"
#include "WLSDetectorConstruction.hh"

#include "WLSActionInitialization.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

// argc holds the number of arguments (including the name) on the command line
// -> it is ONE when only the name is  given !!!
// argv[0] is always the name of the program
// argv[1] points to the first argument, and so on

int main(int argc,char** argv) 
{
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(4);
#else
  // User Verbose output class
  //
  G4VSteppingVerbose::SetInstance(new WLSSteppingVerbose());

  G4int seed = 123;
  if (argc  > 2) seed = atoi(argv[argc-1]);

  // Choose the Random engine and set the seed

  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4Random::setTheSeed(seed);

  G4RunManager * runManager = new G4RunManager;
#endif

  G4String physName = "QGSP_BERT_EMV";

#ifndef WIN32
  G4int c = 0;
  while ((c=getopt(argc,argv,"p")) != -1)
  {
     switch (c)
     {
       case 'p':
         physName = optarg;
         G4cout << "Physics List used is " <<  physName << G4endl;
         break;
       case ':':       /* -p without operand */
         fprintf(stderr, 
                         "Option -%c requires an operand\n", optopt);
         break;
       case '?':
         fprintf(stderr,
                         "Unrecognised option: -%c\n", optopt);
     }
  }
#endif

  // Set mandatory initialization classes

  WLSDetectorConstruction* detector = new WLSDetectorConstruction();

  runManager->SetUserInitialization(detector);

  runManager->SetUserInitialization(new WLSPhysicsList(physName));

  // Set mandatory user action class

  runManager->SetUserInitialization(new WLSActionInitialization(detector));

#ifdef G4VIS_USE
  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager

  G4UImanager * UImanager = G4UImanager::GetUIpointer();

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
  if (argc!=1)
  {
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UImanager->ApplyCommand(command+fileName);
  }
#endif
  else  {
     // Define (G)UI terminal for interactive mode
#ifdef G4UI_USE
     G4UIExecutive * ui = new G4UIExecutive(argc,argv);
#ifdef G4VIS_USE
     UImanager->ApplyCommand("/control/execute init.in");
#endif
     ui->SessionStart();
     delete ui;
#endif
  }

  // job termination

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}
