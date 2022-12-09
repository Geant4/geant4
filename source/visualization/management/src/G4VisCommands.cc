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

// /vis/ top level commands - John Allison  5th February 2001

#include "G4VisCommands.hh"

#include "G4VisManager.hh"
#include "G4UImanager.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4Run.hh"
#include "G4UIsession.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryPoint.hh"
#include "G4RichTrajectory.hh"
#include "G4RichTrajectoryPoint.hh"
#include "G4SmoothTrajectory.hh"
#include "G4SmoothTrajectoryPoint.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4AttDef.hh"

#define G4warn G4cout

////////////// /vis/abortReviewKeptEvents /////////////////////////////

G4VisCommandAbortReviewKeptEvents::G4VisCommandAbortReviewKeptEvents () {
  G4bool omitable;

  fpCommand = new G4UIcmdWithABool("/vis/abortReviewKeptEvents", this);
  fpCommand -> SetGuidance("Abort review of kept events.");
  fpCommand -> SetParameterName("abort", omitable=true);
  fpCommand -> SetDefaultValue(true);
}

G4VisCommandAbortReviewKeptEvents::~G4VisCommandAbortReviewKeptEvents () {
  delete fpCommand;
}

G4String G4VisCommandAbortReviewKeptEvents::GetCurrentValue (G4UIcommand*) {
  return G4String();
}

void G4VisCommandAbortReviewKeptEvents::SetNewValue (G4UIcommand*,
                                                     G4String newValue) {
  fpVisManager->SetAbortReviewKeptEvents(G4UIcommand::ConvertToBool(newValue));
  G4warn << "Type \"continue\" to complete the abort." << G4endl;
}

////////////// /vis/abortReviewPlots /////////////////////////////

G4VisCommandAbortReviewPlots::G4VisCommandAbortReviewPlots () {
  G4bool omitable;

  fpCommand = new G4UIcmdWithABool("/vis/abortReviewPlots", this);
  fpCommand -> SetGuidance("Abort review of plots.");
  fpCommand -> SetParameterName("abort", omitable=true);
  fpCommand -> SetDefaultValue(true);
}

G4VisCommandAbortReviewPlots::~G4VisCommandAbortReviewPlots () {
  delete fpCommand;
}

G4String G4VisCommandAbortReviewPlots::GetCurrentValue (G4UIcommand*) {
  return G4String();
}

void G4VisCommandAbortReviewPlots::SetNewValue (G4UIcommand*,
                                                     G4String newValue) {
  fpVisManager->SetAbortReviewPlots(G4UIcommand::ConvertToBool(newValue));
  G4warn << "Type \"continue\" to complete the abort." << G4endl;
}

////////////// /vis/drawOnlyToBeKeptEvents /////////////////////////////

G4VisCommandDrawOnlyToBeKeptEvents::G4VisCommandDrawOnlyToBeKeptEvents ()
{
  G4bool omitable;

  fpCommand = new G4UIcmdWithABool("/vis/drawOnlyToBeKeptEvents", this);
  fpCommand -> SetGuidance
  ("DURING A RUN draw only those events that have been \"to be kept\" by the user"
   "\nwith G4EventManager::GetEventManager()->KeepTheCurrentEvent() or"
   "\nwith \"/event/keepCurrentEvent\".");
  fpCommand -> SetGuidance(
  "To draw selected events the user should set this flag, then in a user action:"
  "\n  if ( some criterion ) {"
  "\n    G4EventManager::GetEventManager()->KeepTheCurrentEvent();"
  "\n  }");
  fpCommand -> SetParameterName("draw", omitable=true);
  fpCommand -> SetDefaultValue(true);
}

G4VisCommandDrawOnlyToBeKeptEvents::~G4VisCommandDrawOnlyToBeKeptEvents () {
  delete fpCommand;
}

G4String G4VisCommandDrawOnlyToBeKeptEvents::GetCurrentValue (G4UIcommand*) {
  return G4String();
}

void G4VisCommandDrawOnlyToBeKeptEvents::SetNewValue (G4UIcommand*,
                                                  G4String newValue) {
  fpVisManager->SetDrawEventOnlyIfToBeKept(G4UIcommand::ConvertToBool(newValue));
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  if (verbosity < G4VisManager::warnings) {
    if (fpVisManager->GetDrawEventOnlyIfToBeKept()) {
      G4warn << "Only events that have been kept will be drawn." << G4endl;
    } else {
      G4warn << "All events will be drawn." << G4endl;
    }
  }
}

////////////// /vis/enable ///////////////////////////////////////

G4VisCommandEnable::G4VisCommandEnable () {
  G4bool omitable;

  fpCommand = new G4UIcmdWithABool("/vis/enable", this);
  fpCommand -> SetGuidance("Enables/disables visualization system.");
  fpCommand -> SetParameterName("enabled", omitable=true);
  fpCommand -> SetDefaultValue(true);

  fpCommand1 = new G4UIcmdWithoutParameter("/vis/disable", this);
  fpCommand1 -> SetGuidance("Disables visualization system.");
}

G4VisCommandEnable::~G4VisCommandEnable () {
  delete fpCommand;
  delete fpCommand1;
}

G4String G4VisCommandEnable::GetCurrentValue (G4UIcommand*) {
  return G4String();
}

void G4VisCommandEnable::SetNewValue (G4UIcommand* command,
				      G4String newValue) {
  if (command == fpCommand) {
    G4bool enable = G4UIcommand::ConvertToBool(newValue);
    if (enable) fpVisManager->Enable();  // Printing is in vis manager.
    else fpVisManager->Disable();        // Printing is in vis manager.
  } else fpVisManager->Disable();        // Printing is in vis manager.
  // Note: Printing is in vis manager.
}

////////////// /vis/initialize ///////////////////////////////////////

G4VisCommandInitialize::G4VisCommandInitialize ()
{
  fpCommand = new G4UIcmdWithoutParameter("/vis/initialize", this);
  fpCommand -> SetGuidance("Initialise visualisation manager.");
}

G4VisCommandInitialize::~G4VisCommandInitialize () {
  delete fpCommand;
}

void G4VisCommandInitialize::SetNewValue (G4UIcommand*,
					  G4String) {
  fpVisManager->Initialize();
}

////////////// /vis/list ///////////////////////////////////////

G4VisCommandList::G4VisCommandList ()
{
  G4bool omitable;

  fpCommand = new G4UIcmdWithAString("/vis/list", this);
  fpCommand -> SetGuidance("Lists visualization parameters.");
  fpCommand -> SetParameterName("verbosity", omitable=true);
  fpCommand -> SetDefaultValue("warnings");
}

G4VisCommandList::~G4VisCommandList ()
{
  delete fpCommand;
}

G4String G4VisCommandList::GetCurrentValue (G4UIcommand*)
{
  return "";
}

void G4VisCommandList::SetNewValue (G4UIcommand*, G4String newValue)
{
  G4String& verbosityString = newValue;
  G4VisManager::Verbosity verbosity =
    fpVisManager->GetVerbosityValue(verbosityString);

  fpVisManager->PrintAvailableGraphicsSystems(verbosity);
  G4cout << G4endl;
  fpVisManager->PrintAvailableModels(verbosity);
  G4cout << G4endl;
  fpVisManager->PrintAvailableUserVisActions(verbosity);
  G4cout << G4endl;
  fpVisManager->PrintAvailableColours(verbosity);
  G4cout << G4endl;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/vis/scene/list ! " + newValue);
  UImanager->ApplyCommand("/vis/viewer/list ! " + newValue);

  G4cout <<
  "\nAttributes available for modeling and filtering with"
  "\n  \"/vis/modeling/trajectories/create/drawByAttribute\" and"
  "\n  \"/vis/filtering/trajectories/create/attributeFilter\" commands"
  "\nand by picking:"
  << G4endl;
  G4cout
  << *G4TrajectoriesModel().GetAttDefs();
  G4cout
  << *G4RichTrajectory().GetAttDefs()
  << *G4RichTrajectoryPoint().GetAttDefs();
  G4cout
  << *G4SmoothTrajectory().GetAttDefs()
  << *G4SmoothTrajectoryPoint().GetAttDefs();
  G4cout
  << *G4Trajectory().GetAttDefs()
  << *G4TrajectoryPoint().GetAttDefs();

  G4cout <<
  "\nGeometry attributes available for touchables by picking:\n";
  G4cout
  << *G4PhysicalVolumeModel().GetAttDefs();

  if (verbosity < G4VisManager::parameters)
    G4cout <<
  "\nTo get more information, \"/vis/list all\" or use individual commands"
  "\n  such as (use \"ls\" or \"help\"):"
  "\n    /vis/scene/list all all"
  "\n    /vis/viewer/list all all"
  "\n    /vis/modeling/trajectories/list"
  "\n    /vis/filtering/trajectories/list"
	   << G4endl;
}

////////////// /vis/reviewKeptEvents ///////////////////////////////////////

G4VisCommandReviewKeptEvents::G4VisCommandReviewKeptEvents ()
{
  G4bool omitable;

  fpCommand = new G4UIcmdWithAString("/vis/reviewKeptEvents", this);
  fpCommand -> SetGuidance("Review kept events.");
  fpCommand -> SetGuidance
    ("If a macro file is specified, it is executed for each event.");
  fpCommand -> SetGuidance(
  "If a macro file is not specified, each event is drawn to the current"
  "\nviewer.  After each event, the session is paused.  The user may issue"
  "\nany allowed command.  Then enter \"cont[inue]\" to continue to the next"
  "\nevent."
  "\nUseful commands might be:"
  "\n  \"/vis/viewer/...\" to change the view (zoom, set/viewpoint,...)."
  "\n  \"/vis/ogl/export\" to get hard copy."
  "\n  \"/vis/open\" to get alternative viewer."
  "\n  \"/vis/abortReviewKeptEvents\", then \"cont[inue]\", to abort.");
  fpCommand -> SetParameterName("macro-file-name", omitable=true);
  fpCommand -> SetDefaultValue("");
}

G4VisCommandReviewKeptEvents::~G4VisCommandReviewKeptEvents ()
{
  delete fpCommand;
}

G4String G4VisCommandReviewKeptEvents::GetCurrentValue (G4UIcommand*)
{
  return "";
}

void G4VisCommandReviewKeptEvents::SetNewValue (G4UIcommand*, G4String newValue)
{
  if (fpVisManager->GetReviewingKeptEvents()) {
    G4warn <<
      "\"/vis/reviewKeptEvents\" not allowed within an already started review."
      "\n  No action taken."
	   << G4endl;
    return;
  }

  G4String& macroFileName = newValue;
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4RunManager* runManager = G4RunManagerFactory::GetMasterRunManager();
  const G4Run* run         = runManager ? runManager->GetCurrentRun() : nullptr;
  const std::vector<const G4Event*>* events =
    run ? run->GetEventVector() : nullptr;
  size_t nKeptEvents = events ? events->size() : 0;

  if (!nKeptEvents) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<
	"ERROR: G4VisCommandReviewKeptEvents::SetNewValue: No kept events,"
	"\n  or kept events not accessible."
	     << G4endl;
    }
    return;
  }

  G4VViewer* viewer = fpVisManager->GetCurrentViewer();
  if (!viewer) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<
  "ERROR: No current viewer - \"/vis/viewer/list\" to see possibilities."
             << G4endl;
    }
    return;
  }

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn << "ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4int keepControlVerbose = UImanager->GetVerboseLevel();
  G4int newVerbose(0);
  if (keepControlVerbose >= 2 || verbosity >= G4VisManager::confirmations)
    newVerbose = 2;
  UImanager->SetVerboseLevel(newVerbose);
  
  G4VVisManager* keepConcreteInstance = fpVisManager->GetConcreteInstance();
  fpVisManager->Enable();

  // Start on clean view
  UImanager->ApplyCommand("/vis/viewer/rebuild");

  // Event by event refreshing...
  fpVisManager->SetReviewingKeptEvents(true);
  G4bool currentRefreshAtEndOfEvent = pScene->GetRefreshAtEndOfEvent();
  pScene->SetRefreshAtEndOfEvent(true);
  if (macroFileName.empty()) {

    // Draw to viewer and pause session...
    G4UIsession* session = UImanager->GetSession();
    for (size_t i = 0; i < nKeptEvents; ++i) {
      const G4Event* event = (*events)[i];
      if (verbosity >= G4VisManager::warnings) {
	G4warn << "Drawing event : " << event->GetEventID() <<
	  ".  At EndOfEvent, enter any command, then \"cont[inue]\"..."
	       << G4endl;
	static G4bool first = true;
	if (first) {
	  first = false;
	  G4warn <<
  "  Useful commands might be:"
  "\n    \"/vis/scene/add/trajectories\" if not already added."
  "\n    \"/vis/viewer/...\" to change the view (zoom, set/viewpoint,...)."
  "\n    \"/vis/ogl/export\" to get hard copy."
  "\n    \"/vis/open\" to get alternative viewer."
  "\n    \"/vis/abortReviewKeptEvents\", then \"cont[inue]\", to abort."
		 << G4endl;
	}
      }
      fpVisManager->SetRequestedEvent(event);
      UImanager->ApplyCommand("/vis/viewer/rebuild");
      /* The above command forces a rebuild of the scene, including
	 the detector.  This is fine for "immediate" viewers - a
	 refresh requires a rebuild anyway.  But for "stored mode"
	 viewers, you could, in principle, avoid a rebuild of the
	 detector with something like the following:
      sceneHandler->ClearTransientStore();
      viewer->DrawView();
      sceneHandler->DrawEvent(event);
         but this causes mayhem for "immediate" viewers because
         ClearTransientStore issues a DrawView and some curious sort
         of recursion takes place.  For "stored" viewers, the event
         gets drawn but not the eventID, so something odd is happening
         there too.  This needs further investigation - enhanced
         features or a complete re-think.
      */
      UImanager->ApplyCommand("/vis/viewer/flush");
      session->PauseSessionStart("EndOfEvent");
      fpVisManager->SetRequestedEvent(0);
      if (fpVisManager->GetAbortReviewKeptEvents()) break;
    }
    fpVisManager->SetAbortReviewKeptEvents(false);

  } else {

    // Execute macro file...
    for (size_t i = 0; i < nKeptEvents; ++i) {
      const G4Event* event = (*events)[i];
      if (verbosity >= G4VisManager::warnings) {
	G4warn << "Drawing event : " << event->GetEventID()
	       << " with macro file \"" << macroFileName << G4endl;
      }
      fpVisManager->SetRequestedEvent(event);
      UImanager->ApplyCommand("/control/execute " + macroFileName);
      fpVisManager->SetRequestedEvent(0);
    }
  }
  pScene->SetRefreshAtEndOfEvent(currentRefreshAtEndOfEvent);
  fpVisManager->SetReviewingKeptEvents(false);

  if (keepConcreteInstance) fpVisManager->Enable();
  else fpVisManager->Disable();
  UImanager->SetVerboseLevel(keepControlVerbose);
}

////////////// /vis/reviewPlots ///////////////////////////////////////

G4VisCommandReviewPlots::G4VisCommandReviewPlots ()
{
  fpCommand = new G4UIcmdWithoutParameter("/vis/reviewPlots", this);
  fpCommand -> SetGuidance("Review plots.");
  fpCommand -> SetGuidance
  ("Each plot is drawn, one by one, to the current viewer.  After each"
   "\nplot the session is paused.  The user may issue any allowed command."
   "\nThen enter \"cont[inue]\" to continue to the next plot."
   "\nUseful commands might be:"
   "\n  \"/vis/tsg/export\" to get hard copy."
   "\n  \"/vis/abortReviewPlots\", then \"cont[inue]\", to abort.");
}

G4VisCommandReviewPlots::~G4VisCommandReviewPlots ()
{
  delete fpCommand;
}

G4String G4VisCommandReviewPlots::GetCurrentValue (G4UIcommand*)
{
  return "";
}

#include <tools/histo/h1d>
#include <tools/histo/h2d>

namespace {
  template <typename HT>  // tools::histo::h1d, etc
  G4bool ReviewPlots(const G4String& plotType) {  // h1, etc.
    auto visManager = G4VisManager::GetInstance();
    auto ui = G4UImanager::GetUIpointer();
    auto session = ui->GetSession();
    G4bool aborting = false;
    auto keepControlVerbose = ui->GetVerboseLevel();
    ui->SetVerboseLevel(0);
    auto status = ui->ApplyCommand("/analysis/" + plotType + "/getVector");
    ui->SetVerboseLevel(keepControlVerbose);
    if(status==G4UIcommandStatus::fCommandSucceeded) {
      G4String hexString = ui->GetCurrentValues(G4String("/analysis/" + plotType + "/getVector"));
      if(hexString.size()) {
        void* ptr;
        std::istringstream is(hexString);
        is >> ptr;
        auto _v = (const std::vector<HT*>*)ptr;
        auto _n = _v->size();
        for (size_t i = 0; i < _n; ++i) {
          //  Draw then pause session...
          std::ostringstream oss;
          oss << "/vis/plot " << plotType << ' ' << i;
          ui->ApplyCommand(oss.str());
          session->PauseSessionStart("EndOfEvent");
          if (visManager->GetAbortReviewPlots()) {
            aborting = true;
            break;
          }
        }
      }
    }
    return aborting;
  }
}

void G4VisCommandReviewPlots::SetNewValue (G4UIcommand*, G4String)
{
  if (fpVisManager->GetReviewingPlots()) {
    G4warn <<
    "\"/vis/reviewPlots\" not allowed within an already started review."
    "\n  No action taken."
    << G4endl;
    return;
  }

  auto verbosity = fpVisManager->GetVerbosity();

  auto currentViewer = fpVisManager->GetCurrentViewer();
  if (!currentViewer) {
    if (verbosity >= G4VisManager::errors) {
      G4warn <<
      "ERROR: No current viewer - \"/vis/viewer/list\" to see possibilities."
      << G4endl;
    }
    return;
  }

  if (currentViewer->GetName().find("TOOLSSG") == std::string::npos) {
    G4warn <<
    "WARNING: Current viewer not able to draw plots."
    "\n  Try \"/vis/open TSG\", then \"/vis/reviewPlots\" again."
    << G4endl;
    return;
  }

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4warn << "ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  auto ui = G4UImanager::GetUIpointer();

  auto keepControlVerbose = ui->GetVerboseLevel();
  ui->SetVerboseLevel(0);
  auto keepVisVerbose = fpVisManager->GetVerbosity();
  fpVisManager->SetVerboseLevel(G4VisManager::errors);
  auto keepEnable = fpVisManager->IsEnabled();
  fpVisManager->Enable();
  fpVisManager->SetReviewingPlots(true);

  if (ReviewPlots<tools::histo::h1d>("h1")) goto finish;  // Aborting?
  if (ReviewPlots<tools::histo::h2d>("h2")) goto finish;  // Aborting?

finish:
  fpVisManager->SetReviewingPlots(false);
  if (!keepEnable) fpVisManager->Disable();
  fpVisManager->SetVerboseLevel(keepVisVerbose);
  ui->SetVerboseLevel(keepControlVerbose);
}

////////////// /vis/verbose ///////////////////////////////////////

G4VisCommandVerbose::G4VisCommandVerbose () {
  G4bool omitable;

  fpCommand = new G4UIcmdWithAString("/vis/verbose", this);
  for (size_t i = 0; i < G4VisManager::VerbosityGuidanceStrings.size(); ++i) {
    fpCommand -> SetGuidance(G4VisManager::VerbosityGuidanceStrings[i]);
  }
  fpCommand -> SetParameterName("verbosity", omitable=true);
  fpCommand -> SetDefaultValue("warnings");
}

G4VisCommandVerbose::~G4VisCommandVerbose () {
  delete fpCommand;
}

G4String G4VisCommandVerbose::GetCurrentValue (G4UIcommand*) {
  return G4String();
}

void G4VisCommandVerbose::SetNewValue (G4UIcommand*,
				       G4String newValue) {
  G4VisManager::Verbosity verbosity =
    fpVisManager->GetVerbosityValue(newValue);
  fpVisManager->SetVerboseLevel(verbosity);
  // Always prints whatever the verbosity...
  G4cout << "Visualization verbosity changed to "
	 << G4VisManager::VerbosityString(verbosity) << G4endl;
}
