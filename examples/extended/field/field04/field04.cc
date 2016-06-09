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
// $Id: field04.cc,v 1.8.2.1 2009/08/11 09:47:40 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02-patch-04 $
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

#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

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

// argc holds the number of arguments (including the name) on the command line
// -> it is ONE when only the name is  given !!!
// argv[0] is always the name of the program
// argv[1] points to the first argument, and so on

int main(int argc,char** argv) 
{
  G4bool useUItcsh = true;

  G4String physicsList = "QGSP_BERT";

  G4int seed = 123;
  if (argc  > 2) seed = atoi(argv[argc-1]);
  // if (seed == 123) seed = time(0);

  // Choose the Random engine

  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(seed);
  CLHEP::HepRandom::showEngineStatus();

#ifndef WIN32
  G4int c = 0;
  while ((c=getopt(argc,argv,"p:t")) != -1)
  {
     switch (c)
     {
       case 'p':
         physicsList = optarg;
         G4cout << "Physics List used is " <<  physicsList << G4endl;
         break;
       case 't': // Don't use a tcsh-style command line interface
         useUItcsh = false;
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
  runManager->SetUserInitialization(new F04PhysicsList(physicsList));
  
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

  G4UImanager * UI = G4UImanager::GetUIpointer();  

#ifndef WIN32
  G4int optmax = argc;
  if (argc > 2)  { optmax = optmax-1; }

  if (optind < optmax)
  {
     G4String command = "/control/execute ";
     for ( ; optind < optmax; optind++)
     {
         G4String macroFilename = argv[optind];
         UI->ApplyCommand(command+macroFilename);
     }
  }
  else
  {
     // Define (G)UI terminal for interactive mode
     G4UIsession * session = 0;
     if (useUItcsh)
     {
        // G4UIterminal is a terminal with tcsh-like control.
        session = new G4UIterminal(new G4UItcsh);
     }
     else
     {
        // G4UIterminal is a (dumb) terminal.
        session = new G4UIterminal();
     }
     session->SessionStart();
     delete session;
  }
#else  // Simple UI for Windows runs, no possibility of additional arguments
  if (argc==1)   // Define UI terminal for interactive mode  
  { 
     G4UIsession * session = new G4UIterminal;
     session->SessionStart();
     delete session;
  }
  else           // Batch mode
  { 
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UI->ApplyCommand(command+fileName);
  }
#endif      
  // job termination

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}
