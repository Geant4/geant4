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
// $Id: field04.cc,v 1.12 2010-10-25 00:49:58 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// --------------------------------------------------------------
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

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "F04PhysicsList.hh"
#include "F04DetectorConstruction.hh"
#include "F04PrimaryGeneratorAction.hh"

#include "F04RunAction.hh"
#include "F04EventAction.hh"
#include "F04TrackingAction.hh"
#include "F04SteppingAction.hh"
#include "F04StackingAction.hh"
#include "F04SteppingVerbose.hh"

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
  G4String physName = "QGSP_BERT_EMV";

  G4int seed = 123;
  if (argc  > 2) seed = atoi(argv[argc-1]);

  // Choose the Random engine

  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(seed);

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

  // My Verbose output class

  G4VSteppingVerbose::SetInstance(new F04SteppingVerbose);
  
  // Construct the default run manager

  G4RunManager * runManager = new G4RunManager;

  // Set mandatory initialization classes

  F04DetectorConstruction* detector = new F04DetectorConstruction();

  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new F04PhysicsList(physName));
  
#ifdef G4VIS_USE

  // visualization manager

  G4VisManager* visManager = new G4VisExecutive();
  visManager->Initialize();

#endif

  // Set mandatory user action class 

  runManager->SetUserAction( new F04PrimaryGeneratorAction(detector) );

  F04RunAction* runAction = new F04RunAction();
  F04EventAction* eventAction = new F04EventAction(runAction);

  runManager->SetUserAction(runAction);
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction( new F04TrackingAction() );
  runManager->SetUserAction( new F04SteppingAction() );
  runManager->SetUserAction( new F04StackingAction() );

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
     UImanager->ApplyCommand("/control/execute vis.mac");     
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
