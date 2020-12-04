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
/// \file parallel/ThreadsafeScorers/ts_scorers.cc
/// \brief Main of the ThreadsafeScorers example
//
//
//
//
/// ts_scorers example shows how to use global scorers. The benefit of using
///     global scorers in memory-savings for problems with very large amounts
///     of scoring volumes. Additionally, the global scorers are more precise
///     w.r.t. the serial solution because of the lack of compounding
///     round-off error from multiple threads
///
/// In this example, the global scorers are implemented as static member
///     variables in TSRun because TSRun is thread-local. The G4atomic
///     class is the core of the thread-safe scorers and can be uses
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Types.hh"

#include "G4RunManagerFactory.hh"

#include "G4Threading.hh"

#include "Randomize.hh"

// User Defined Classes
#include "TSActionInitialization.hh"
#include "TSDetectorConstruction.hh"
#include "TSPhysicsList.hh"

#include "G4TiMemory.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4Track.hh"
#include "G4Step.hh"

// for std::system(const char*)
#include <cstdlib>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void message(G4RunManager* runmanager)
{
  G4MTRunManager* man = dynamic_cast<G4MTRunManager*>(runmanager);
  if(man)
  {
    man->SetNumberOfThreads(G4Threading::G4GetNumberOfCores());
    G4cout << "\n\n\t--> Running in multithreaded mode with "
           << man->GetNumberOfThreads() << " threads\n\n"
           << G4endl;
  }
  else
  {
    G4cout << "\n\n\t--> Running in serial mode\n\n" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  // initialize timemory
  G4Profiler::Configure(argc, argv);

  G4String macro;
  if(argc > 1)
    macro = argv[argc - 1];

  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = 0;
  if(macro.empty())
    ui = new G4UIExecutive(argc, argv);

  // Set the random seed
  CLHEP::HepRandom::setTheSeed(1245214UL);

#if defined(GEANT4_USE_TIMEMORY)
  // The following exists for:
  // - G4ProfileType::Run
  // - G4ProfileType::Event
  // - G4ProfileType::Track
  // - G4ProfileType::Step
  // - G4ProfileType::User
  //
  using TrackProfilerConfig = G4ProfilerConfig<G4ProfileType::Track>;
  using TrackTool           = typename TrackProfilerConfig::type;

  TrackProfilerConfig::GetQueryFunctor() = [](const G4Track* _track) {
    // only profile if _track != nullptr and dynamic-profiler != nullptr
    // and /profiler/track/enable is true
    //
    return G4Profiler::GetEnabled(G4ProfileType::Track) && _track &&
           _track->GetDynamicParticle();
  };

  TrackProfilerConfig::GetLabelFunctor() = [](const G4Track* _track) {
    // create a label for the profiling entry. This can be customized
    // to include and information necessary in the returning string
    auto pdef = _track->GetDynamicParticle()->GetParticleDefinition();
    static std::string _prefix = "G4Track/";
    return _prefix + pdef->GetParticleName();
  };

  // env option to display track profiles as a hierarchy
  bool track_tree = tim::get_env<bool>("G4PROFILER_TRACK_TREE", true);
  // env option to enable timeline entries (every entry is unique, HUGE amount
  // of data!)
  bool track_time = tim::get_env<bool>("G4PROFILER_TRACK_TIMELINE", false);
  // default scope is tree
  auto _scope = tim::scope::config{};
  if(track_tree == false)
    _scope += tim::scope::flat{};
  if(track_time == true)
    _scope += tim::scope::timeline{};
  TrackProfilerConfig::GetToolFunctor() = [=](const std::string& _label) {
    // Configure the profiling tool for a given label. By default,
    // G4Track and G4Step tools are "flat profiles" but this can be disabled
    // to include tree
    return new TrackTool(_label, _scope);
  };
#endif

  G4RunManager* runmanager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Tasking);

  message(runmanager);

  runmanager->SetUserInitialization(new TSDetectorConstruction);

  runmanager->SetUserInitialization(new TSPhysicsList);

  runmanager->SetUserInitialization(new TSActionInitialization);

  runmanager->Initialize();

  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if(!ui)
  {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + macro);
  }
  else
  {
    ui->SessionStart();
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !
  delete visManager;
  delete runmanager;

  return 0;
}
