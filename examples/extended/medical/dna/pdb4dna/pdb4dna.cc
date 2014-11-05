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
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $Id$
//
/// \file pdb4dna.cc
/// \brief Main program of the pdb4dna example

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

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

#include "CommandLineParser.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace G4DNAPARSER;
CommandLineParser* parser(0);

void Parse(int& argc, char** argv);


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  //////////
  // Parse options given in commandLine
  //
  Parse(argc, argv);

  // Set the Seed
  CLHEP::RanecuEngine defaultEngine(1234567);
  G4Random::setTheEngine(&defaultEngine);

  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  //////////
  // Construct the run manager according to whether MT is activated or not
  //
  G4RunManager* runManager(0);
  Command* commandLine(0);

  if ((commandLine = parser->GetCommandIfActive("-mt")))
  {
#ifdef G4MULTITHREADED

    runManager= new G4MTRunManager;

    if(commandLine->GetOption().empty())
    {
      ((G4MTRunManager*)runManager)->SetNumberOfThreads(1);
    }
    else
    {
      int nThreads = G4UIcommand::ConvertToInt(commandLine->GetOption());
      ((G4MTRunManager*)runManager)->SetNumberOfThreads(nThreads);
    }
#else
    G4cout << "WARNING : the -mt command line option has not effect since you "
        "seam to have compile Geant4 without the G4MULTITHREADED flag"
           << G4endl;
    runManager = new G4RunManager();

#endif
  }
  else
  {
    runManager = new G4RunManager;
  }
  
  // Set user classes
  //
  runManager->SetUserInitialization(new DetectorConstruction());
  runManager->SetUserInitialization(new PhysicsList);
  runManager->SetUserInitialization(new ActionInitialization());

  // Initialize G4 keRnel
  //
  // runManager->Initialize();

#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ((commandLine = parser->GetCommandIfActive("-mac")) ) {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+commandLine->GetOption());
  }
  else
  {
    UImanager->ApplyCommand("/control/execute init.mac");
 
    // interactive mode : define UI session
    if ((commandLine = parser->GetCommandIfActive("-gui")))
    {
      G4UIExecutive* ui = new G4UIExecutive(argc, argv,
                                            commandLine->GetOption());

      if(parser->GetCommandIfActive("-novis") == 0) 
      // visualization is used by default
      {
        if ((commandLine = parser->GetCommandIfActive("-vis")))
        // select a visualization driver if needed (e.g. HepFile)
        {
          UImanager->ApplyCommand(G4String("/vis/open ")+commandLine->GetOption()); 
        } 
        else
        // by default OGL is used
        {
          UImanager->ApplyCommand("/vis/open OGL 600x600-0+0");
        }
        UImanager->ApplyCommand("/control/execute vis.mac");
      }
      if (ui->IsGUI())
        UImanager->ApplyCommand("/control/execute gui.mac");
      ui->SessionStart();
      delete ui;
    }
   else
   {
     if ((commandLine = parser->GetCommandIfActive("-vis")))
     {
        UImanager->ApplyCommand(G4String("/vis/open ")+commandLine->GetOption());
        UImanager->ApplyCommand("/control/execute vis.mac"); 
     }
   }
  }

  // Job termination

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void Parse(int& argc, char** argv)
{
  //////////
  // Parse options given in commandLine
  //
  parser = CommandLineParser::GetParser();
#ifdef G4UI_USE
  parser->AddCommand(
      "-gui", Command::OptionNotCompulsory,
      "Select geant4 UI or just launch a geant4 terminal session", "qt");
#endif
  parser->AddCommand("-mac", Command::WithOption, "Give a mac file to execute",
                     "pdb4dna.in");
// You cann your own command, as for instance:
//  parser->AddCommand("-seed", Command::WithOption,
//                     "Give a seed value in argument to be tested", "seed");
// it is then up to you to manage this option
  parser->AddCommand("-mt", Command::OptionNotCompulsory,
                     "Launch in MT mode (events computed in parallel)",
                     "2");

  parser->AddCommand("-vis", Command::WithOption, "Select a visualization driver",
                     "OGL 600x600-0+0");


  parser->AddCommand("-novis", Command::WithoutOption, "Deactivate visualization when using GUI");


  //////////
  // If -h or --help is given in option : print help and exit
  //
  if (parser->Parse(argc, argv) != 0) // help is being printed
  {
    // if you are using ROOT, create a TApplication in this condition in order
    // to print the help from ROOT as well
    CommandLineParser::DeleteInstance();
    std::exit(0);
  }

  ///////////
  // Kill application if wrong argument in command line
  //
  if (parser->CheckIfNotHandledOptionsExists(argc, argv))
  {
    // if you are using ROOT, you should initialise your TApplication
    // before this condition
    abort();
  }
}

