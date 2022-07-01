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
// G4ProfilerMessenger implementation
//
// Author: J.Madsen, 12 November 2020
// --------------------------------------------------------------------

#include "G4ProfilerMessenger.hh"

#include "G4Profiler.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4TiMemory.hh"

// --------------------------------------------------------------------

using Type = G4ProfileType;

G4ProfilerMessenger::G4ProfilerMessenger()
{
  profileDirectory = new G4UIdirectory("/profiler/");
  profileDirectory->SetGuidance("Profiler controls.");

  profileOutputDirectory = new G4UIdirectory("/profiler/output/");
  profileOutputDirectory->SetGuidance(
    "Control the output modes of the profiler.");

#define CREATE_DIR(IDX, DIR, GUIDANCE)                                         \
  profileTypeDirs.at(IDX) = new G4UIdirectory(DIR);                            \
  profileTypeDirs.at(IDX)->SetGuidance(GUIDANCE)

  CREATE_DIR(Type::Run, "/profiler/run/",
             "Profiler controls at the G4Run level");
  CREATE_DIR(Type::Event, "/profiler/event/",
             "Profiler controls at the G4Event level");
  CREATE_DIR(Type::Track, "/profiler/track/",
             "Profiler controls at the G4Track level");
  CREATE_DIR(Type::Step, "/profiler/step/",
             "Profiler controls at the G4Step level");
  CREATE_DIR(Type::User, "/profiler/user/",
             "Profiler controls within user code");

#define SET_ENABLED_CMD(IDX, CMD, CMDLINE, DEFAULT, GUIDANCE)                  \
  profileEnableCmds.at(IDX).second = CMDLINE;                                  \
  profileEnableCmds.at(IDX).first  = new G4UIcmdWithABool(CMD, this);          \
  profileEnableCmds.at(IDX).first->SetDefaultValue(DEFAULT);                   \
  profileEnableCmds.at(IDX).first->SetGuidance(GUIDANCE);                      \
  profileEnableCmds.at(IDX).first->AvailableForStates(G4State_PreInit,         \
                                                      G4State_Idle)

  SET_ENABLED_CMD(Type::Run, "/profiler/run/enable", "run", true,
                  "Record metrics for each G4Run");
  SET_ENABLED_CMD(Type::Event, "/profiler/event/enable", "event", true,
                  "Record metrics for each G4Event");
  SET_ENABLED_CMD(Type::Track, "/profiler/track/enable", "track", false,
                  "Record metrics for each G4Track");
  SET_ENABLED_CMD(Type::Step, "/profiler/step/enable", "step", false,
                  "Record metrics for each G4Step");
  SET_ENABLED_CMD(Type::User, "/profiler/user/enable", "user", true,
                  "Record metrics for user specified profiling instances");

#define SET_COMPONENTS_CMD(IDX, CMD, CMDLINE, DEFAULTS, GUIDANCE)              \
  profileCompCmds.at(IDX).second = CMDLINE;                                    \
  profileCompCmds.at(IDX).first  = new G4UIcmdWithAString(CMD, this);          \
  profileCompCmds.at(IDX).first->SetDefaultValue(DEFAULTS);                    \
  profileCompCmds.at(IDX).first->SetGuidance(GUIDANCE);                        \
  profileCompCmds.at(IDX).first->AvailableForStates(G4State_PreInit,           \
                                                    G4State_Idle)

  G4String comps = "wall_clock, cpu_clock, cpu_util, peak_rss";

  SET_COMPONENTS_CMD(
    Type::Run, "/profiler/run/components", "--run-components", comps,
    "Measurment types to record for each G4Run (see `timemory-avail -s`)");
  SET_COMPONENTS_CMD(
    Type::Event, "/profiler/event/components", "--event-components", comps,
    "Measurment types to record for each G4Event (see `timemory-avail -s`)");
  SET_COMPONENTS_CMD(
    Type::Track, "/profiler/track/components", "--track-components", comps,
    "Measurment types to record for each G4Track (see `timemory-avail -s`)");
  SET_COMPONENTS_CMD(
    Type::Step, "/profiler/step/components", "--step-components", comps,
    "Measurment types to record for each G4Step (see `timemory-avail -s`)");
  SET_COMPONENTS_CMD(Type::User, "/profiler/user/components",
                     "--user-components", comps,
                     "Measurment types to record for user specified profiling "
                     "instances (see `timemory-avail -s`)");

#define SET_OUTPUT_CMD(CMD, CMDLINE, DEFAULT, GUIDANCE)                        \
  profileGeneralCmds.push_back({ new G4UIcmdWithABool(CMD, this), CMDLINE });  \
  profileGeneralCmds.back().first->SetDefaultValue(DEFAULT);                   \
  profileGeneralCmds.back().first->SetGuidance(GUIDANCE);                      \
  profileGeneralCmds.back().first->AvailableForStates(G4State_PreInit,         \
                                                      G4State_Idle)

  SET_OUTPUT_CMD("/profiler/output/dart", "--dart", false,
                 "Enabled Dart output (CTest/CDash data tracking)");
  SET_OUTPUT_CMD("/profiler/output/json", "--json", true,
                 "Enabled JSON output");
  SET_OUTPUT_CMD("/profiler/output/text", "--text", true,
                 "Enabled text output");
  SET_OUTPUT_CMD("/profiler/output/cout", "--cout", false,
                 "Enabled output to console");
  SET_OUTPUT_CMD("/profiler/output/plot", "--plot", false,
                 "Enabled plotting JSON output");

  SET_OUTPUT_CMD("/profiler/tree", "--tree", true,
                 "Display the results as a call-stack hierarchy.");
  SET_OUTPUT_CMD("/profiler/flat", "--flat", false,
                 "Display the results as a flat call-stack");
  SET_OUTPUT_CMD("/profiler/timeline", "--timeline", false,
                 "Do not merge duplicate entries at the same call-stack "
                 "position. May be combined with tree or flat profiles.");
  SET_OUTPUT_CMD(
    "/profiler/per_thread", "--per-thread", false,
    "Display the results for each individual thread (default: aggregation)");
  SET_OUTPUT_CMD(
    "/profiler/per_event", "--per-event", false,
    "Display the results for each individual G4event (default: aggregation)");
}

// --------------------------------------------------------------------

G4ProfilerMessenger::~G4ProfilerMessenger()
{
  delete profileDirectory;
  delete profileOutputDirectory;
  for(auto& itr : profileTypeDirs)
  {
    delete itr;
  }
  for(auto& itr : profileEnableCmds)
  {
    delete itr.first;
  }
  for(auto& itr : profileGeneralCmds)
  {
    delete itr.first;
  }
  for(auto& itr : profileCompCmds)
  {
    delete itr.first;
  }
}

// --------------------------------------------------------------------

void G4ProfilerMessenger::SetNewValue(G4UIcommand* command, G4String value)
{
  for(size_t i = 0; i < static_cast<size_t>(G4ProfileType::TypeEnd); ++i)
  {
    G4UIcmdWithABool* ui = profileEnableCmds.at(i).first;
    if(command == ui)
    {
      G4Profiler::SetEnabled(i, ui->GetNewBoolValue(value));
      return;
    }
  }

  // pass the commands to the timemory argparser
  std::vector<std::string> command_line = { "G4ProfilerMessenger" };

  for(auto& itr : profileGeneralCmds)
  {
    G4UIcmdWithABool* ui = itr.first;
    if(command == ui)
    {
      command_line.push_back(itr.second);
      command_line.push_back(value);
      break;
    }
  }

  for(auto& itr : profileCompCmds)
  {
    G4UIcmdWithAString* ui = itr.first;
    if(command == ui)
    {
      command_line.push_back(itr.second);
#if defined(GEANT4_USE_TIMEMORY)
      for(auto vitr : tim::delimit(value, ", ;"))
        command_line.push_back(vitr);
#endif
      break;
    }
  }

  if(command_line.size() > 1)
  {
    G4Profiler::Configure(command_line);
  }
}

// --------------------------------------------------------------------
