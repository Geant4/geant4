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

#include "LaunchG4.hh"
#include "Parser.hh"
#include "CLHEP/Random/Random.h"
#include "Randomize.hh"
#include "G4UImanager.hh"
#include <string>
#include <time.h>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Parser* parser = 0;
LaunchG4 *g4 = 0; // Must be created before calling parser->Parse(argc,argv);
G4int seed = 0;

void SetSeed();
unsigned TokenizeString(const std::string& i_source,
                        const std::string& i_seperators,
                        bool i_discard_empty_tokens,
                        std::vector<std::string>& o_tokens);

int main(int argc, char* argv[])
{
  parser = Parser::GetParser();
  g4 = new LaunchG4(); // Must be created before calling parser->Parse(argc,argv);

  G4bool g4session = false;
  G4String macFile = "";
  G4bool useCustomMacFile = false;
  G4bool chemistryflag = true;

  int firstArgc = argc;
  parser->AddCommand(
      "-gui", OptionNotCompulsory,
      "Select geant4 UI or just launch a geant4 terminal session", "qt");
  parser->AddCommand("-mac", WithOption, "Give a mac file to execute",
                     "macFile.mac");
  parser->AddCommand("-seed", WithOption,
                     "Give a seed value in argument to be tested", "seed");
  parser->AddCommand(
      "-mt",
      OptionNotCompulsory,
      "Launch in MT mode (events computed in parallel, NOT RECOMMANDED WITH CHEMISTRY)");
  parser->AddCommand("-cluster", WithoutOption,
                     "Specify that you are working on the CENBG cluster");

  if (parser->Parse(argc, argv) != 0) // First initialize LaunchG4 and root then parse options
  {
    Parser::DeleteInstance();
    delete g4;
    std::exit(0);
  }

  // Get the pointer to the User Interface manager
//	G4UImanager* UI = G4UImanager::GetUIpointer();

  if (firstArgc > 1) // Batch mode
  {
    Command* commandLine(0);
    if ((commandLine = parser->GetCommandIfActive("-mac")))
    {
      useCustomMacFile = true;
      macFile = commandLine->fOption;
    }
    if ((commandLine = parser->GetCommandIfActive("-seed")))
    {
      seed = atoi(commandLine->fOption.c_str());
    }
    if ((commandLine = parser->GetCommandIfActive("-gui")))
    {
      g4session = true;
    }

    // remove handled arguments from argument array
    int j = 0;
    for (int i = 0; i < argc; i++)
    {
      if (strcmp(argv[i], ""))
      {
        argv[j] = argv[i];
        j++;
      }
    }
    argc = j;
  }

  SetSeed();
  G4cout << "Seed used : " << seed << G4endl;
  // Choose the Random engine
  // CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  g4->Initialize(chemistryflag);

  if (g4session)
  {
    if (Command* commandLine = parser->GetCommandIfActive("-gui"))
    {
      g4->NewSession(argc, argv, commandLine->fOption);
    }
    else
    {
      g4->NewSession(argc, argv);
    }
  }
  else if (useCustomMacFile == false)
  {
    g4->RunSimu("commands.mac");
    //g4 -> NewSession(argc,argv,"qt");
    g4session = true;
  }

  //____________________________________________
  // Kill application if wrong argument in command line
  if (argc > 0)
  {
    G4bool kill = false;
    for (G4int i = 1; i < argc; i++)
    {
      if (strcmp(argv[i], ""))
      {
        kill = true;
        G4cerr << "Unknown argument : " << argv[i] << "\n";
      }
    }
    if (kill)
    {
      G4cerr << "Le programme " << argv[0]
             << " va être tué par le main : pb d'argument" << G4endl;
      abort(); // KILL APPLICATION
    }
  }

  if (useCustomMacFile)
  {
    //____________________________________
    // Loop over mac files in command line
//        if(macFile.empty() == false)
    {
      G4cout << "macFile = " << macFile << G4endl;

      vector<std::string> tokens;
      TokenizeString(macFile," ,{}",true,tokens);
      vector<std::string>::iterator it = tokens.begin();

      for(; it != tokens.end(); it++)
      {
        g4 -> RunSimu(*it);
      }
    }
    //____________________________________
  }

  if (g4session)
  {
    g4->StartSession();
  }

#ifdef USE_ROOT
  delete root;
#endif

  delete g4;

  delete parser;
  return 0;
}

void SetSeed()
{
  // Choose the Random engine
  // CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  if (seed == 0) // If no seed given in argument, setup the seed
  {
    int jobID_int = 0;

    //____________________________________
    // In case on cluster
    Command* commandLine(0);
    if ((commandLine = parser->GetCommandIfActive("-cluster")))
    {
      const char * env = getenv("PBS_JOBID");

      if (env)
      {
        string buffer(env);
        string jobID_string = buffer.substr(0, buffer.find("."));

        if (jobID_string.find("-") < jobID_string.length())
        {
          jobID_string = jobID_string.substr(jobID_string.find("-") + 1,
                                             jobID_string.length());
        }

        jobID_int = atoi(jobID_string.c_str());
      }
    }
    // end cluster
    //____________________________________
    seed = ((int) time(NULL)) + jobID_int;
    //g4->SetJobID(jobID_int);
    G4Random::setTheSeed(seed);
  }
  else
  {
    G4Random::setTheSeed(seed);
  }
}

unsigned TokenizeString(const std::string& i_source,
                        const std::string& i_seperators,
                        bool i_discard_empty_tokens,
                        std::vector<std::string>& o_tokens)
{
  size_t prev_pos = 0;
  size_t pos = 0;
  unsigned number_of_tokens = 0;
  o_tokens.clear();
  pos = i_source.find_first_of(i_seperators, pos);
  while (pos != std::string::npos)
  {
    std::string token = i_source.substr(prev_pos, pos - prev_pos);
    if (!i_discard_empty_tokens || token != "")
    {
      o_tokens.push_back(i_source.substr(prev_pos, pos - prev_pos));
      number_of_tokens++;
    }

    pos++;
    prev_pos = pos;
    pos = i_source.find_first_of(i_seperators, pos);
  }

  if (prev_pos < i_source.length())
  {
    o_tokens.push_back(i_source.substr(prev_pos));
    number_of_tokens++;
  }

  return number_of_tokens;
}

