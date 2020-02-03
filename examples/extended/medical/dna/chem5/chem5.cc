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
// J. Comput. Phys. 274 (2014) 841-882
// Phys. Med. Biol. 63(10) (2018) 105014-12pp 
// The Geant4-DNA web site is available at http://geant4-dna.org
//
//
/// \file chem5.cc
/// \brief Chem5 example

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4DNAChemistryManager.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

#include "CommandLineParser.hh"

/*
 * WARNING : Geant4 was initially not intended for this kind of application
 * This code is delivered as a prototype
 * We will be happy to hear from you, do not hesitate to send your feedback and
 * communicate on the difficulties you may encounter
 * The user interface may change in the next releases since a reiteration of the
 * code has started
 */

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace G4DNAPARSER;
CommandLineParser* parser(0);
long seed = 0;

unsigned int noise();
void SetSeed();

void Parse(int& argc, char** argv);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  // Parse options given in commandLine
  Parse(argc, argv);
  Command* commandLine(0);
  SetSeed();
  
  // Construct the run manager according to whether MT is activated or not
  //
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager= new G4MTRunManager();
  
  if((commandLine = parser->GetCommandIfActive("-mt")))
  {
    int nThreads = 2;
    if(commandLine->GetOption() == "NMAX")
    {
      nThreads = G4Threading::G4GetNumberOfCores();
    }
    else
    {
      nThreads = G4UIcommand::ConvertToInt(commandLine->GetOption());
    }
    
    runManager->SetNumberOfThreads(nThreads);
  }
  
  G4cout << "**************************************************************"
         << "******\n===== Chem5 is started with "
         << runManager->GetNumberOfThreads() << " of "
         << G4Threading::G4GetNumberOfCores()
         << " available threads =====\n\n*************************************"
         <<"*******************************" 
         << G4endl;
#else
  G4RunManager* runManager = new G4RunManager();
#endif
  
  // Set mandatory initialization classes
  runManager->SetUserInitialization(new PhysicsList());
  runManager->SetUserInitialization(new DetectorConstruction());
  runManager->SetUserInitialization(new ActionInitialization());
  
  // Initialize G4 kernel
  runManager->Initialize();
  
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive();
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4UIExecutive* ui(0);
  
  // interactive mode : define UI session
  if((commandLine = parser->GetCommandIfActive("-gui")))
  {
    ui = new G4UIExecutive(argc, argv, commandLine->GetOption());
    
    if (ui->IsGUI()) UImanager->ApplyCommand("/control/execute gui.mac");
    
    if (parser->GetCommandIfActive("-novis") == 0)
    {
      // visualization is used by default
      if ((commandLine = parser->GetCommandIfActive("-vis")))
      {
        // select a visualization driver if needed (e.g. HepFile)
        UImanager->ApplyCommand
        (G4String("/vis/open ") + commandLine->GetOption());
      }
      else
      {
        // by default OGL is used
        UImanager->ApplyCommand("/vis/open OGL 800x600-0+0");
      }
      UImanager->ApplyCommand("/control/execute vis.mac");
    }
  }
  else
  {
    // to be use visualization file (= store the visualization into
    // an external file:
    // ASCIITree ;  DAWNFILE ; HepRepFile ; VRML(1,2)FILE ; gMocrenFile ...
    if((commandLine = parser->GetCommandIfActive("-vis")))
    {
      UImanager->ApplyCommand(G4String("/vis/open ")
                              + commandLine->GetOption());
      UImanager->ApplyCommand("/control/execute vis.mac");
    }
  }
  
  if((commandLine = parser->GetCommandIfActive("-mac")))
  {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + commandLine->GetOption());
  }
  else
  {
    UImanager->ApplyCommand("/control/execute beam.in");
  }

  if((commandLine = parser->GetCommandIfActive("-gui")))
  {
    ui->SessionStart();
    delete ui;
  }
  
  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !
  delete visManager;
  delete runManager;
  CommandLineParser::DeleteInstance();
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool IsBracket(char c)
{
  switch(c)
  {
    case '[':
    case ']':
      return true;
    default:
      return false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SetSeed()
{  
  Command* commandLine(0);
  
  if((commandLine = parser->GetCommandIfActive("-seed")))
  {
    seed = atoi(commandLine->GetOption().c_str());
  }
  
  if(seed == 0) // If no seed given in argument, setup the seed
  {
    long jobID_int = 0;
    long noice = 0;
    
    //____________________________________
    // In case on cluster
    if((commandLine = parser->GetCommandIfActive("-cluster")))
    {
      noice = labs((long) noise());
      
      const char * env = std::getenv("PBS_JOBID");

      if(env)
      {
        G4String buffer(env);
        G4String jobID_string = buffer.substr(0,  buffer.find("."));
        jobID_string.erase(std::remove_if(jobID_string.begin(),
                                          jobID_string.end(),
                                          &IsBracket),
                                          jobID_string.end());
        jobID_int = atoi(jobID_string.c_str());
      }
      else
      {
        env = std::getenv("SGE_TASK_ID");
        if(env) jobID_int = atoi(env);
      }
    } // end cluster
    
    //____________________________________
    seed = ((long) time(NULL)) + jobID_int + noice;
  }
  
  G4cout << "Seed used : " << seed << G4endl;
  G4Random::setTheEngine(new CLHEP::MixMaxRng());
  G4Random::setTheSeed(seed);

  // Choose the Random engine
  // CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

unsigned int noise()
{
#if defined(WIN32)|| defined(_WIN32)|| defined(__WIN32)&&!defined(__CYGWIN__)
  // TODO: MS Win
  return std::time(0);
#else
  unsigned int random_seed, random_seed_a, random_seed_b;
  std::ifstream file ("/dev/urandom", std::ios::binary);
  if (file.is_open())
  {
    char * memblock;
    int size = sizeof(int);
    memblock = new char [size];
    file.read (memblock, size);
    file.close();
    random_seed_a = *reinterpret_cast<int*>(memblock);
    delete[] memblock;
  }// end if
  else
  {
    random_seed_a = 0;
  }
  random_seed_b = std::time(0);
  random_seed = random_seed_a xor random_seed_b;
  return random_seed;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Parse(int& argc, char** argv)
{
  //////////
  // Parse options given in commandLine
  //
  parser = CommandLineParser::GetParser();
  
  parser->AddCommand("-gui", Command::OptionNotCompulsory,
                     "Select geant4 UI or just launch a geant4 terminal"
                     "session",
                     "qt");
                     
  parser->AddCommand("-mac", Command::WithOption, "Give a mac file to execute",
                     "macFile.mac");

  parser->AddCommand("-seed",
                     Command::WithOption,
                     "Give a seed value in argument to be tested", "seed");
  
#ifdef G4MULTITHREADED
  parser->AddCommand("-mt",Command::WithOption,
                     "Launch in MT mode (events computed in parallel,"
                     " NOT RECOMMANDED WITH CHEMISTRY)", "2");
#endif
  
  parser->AddCommand("-chemOFF", Command::WithoutOption,
                     "Deactivate chemistry");
                     
  parser->AddCommand("-vis", Command::WithOption,
                     "Select a visualization driver", "OGL 600x600-0+0");
                     
  parser->AddCommand("-novis", Command::WithoutOption,
                     "Deactivate visualization when using GUI");
  
  parser->AddCommand("-cluster", Command::WithoutOption,
                     "Launch the code on a cluster, avoid dupplicated seeds");
  
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
  if(parser->CheckIfNotHandledOptionsExists(argc, argv))
  {
    // if you are using ROOT, you should initialise your TApplication
    // before this condition
    std::exit(0);
  }
}
