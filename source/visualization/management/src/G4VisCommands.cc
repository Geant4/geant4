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
// $Id: G4VisCommands.cc,v 1.20 2006/11/26 15:49:10 allison Exp $
// GEANT4 tag $Name: geant4-08-02 $

// /vis/ top level commands - John Allison  5th February 2001

#include "G4VisCommands.hh"

#include "G4VisManager.hh"
#include "G4UImanager.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4UIsession.hh"

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

  fpVisManager->PrintAvailableGraphicsSystems();
  G4cout << G4endl;
  fpVisManager->PrintAvailableModels(verbosity);
  G4cout << G4endl;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand(G4String("/vis/viewer/list ! ") + verbosityString);
  if (verbosity < G4VisManager::parameters)
    G4cout <<
  "\nTo get more information, \"/vis/list all all\" or use individual commands"
  "\n  such as (use \"ls\" or \"help\"):"
  "\n    /vis/viewer/list"
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
  fpCommand -> SetGuidance
    ("If a macro file is not specified, each event is drawn to the current"
     "\nviewer.  After each event, the session is paused.  The user may issue"
     "\nany allowed command.  Then enter \"continue\" to continue to the next"
     "\nevent.");
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
  G4String& macroFileName = newValue;
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4RunManager* runManager = G4RunManager::GetRunManager();
  const G4Run* run = runManager? runManager->GetCurrentRun(): 0;
  const std::vector<const G4Event*>* events = run? run->GetEventVector(): 0;
  size_t nKeptEvents = events? events->size(): 0;

  if (!nKeptEvents) {
    if (verbosity >= G4VisManager::errors) {
      G4cout <<
	"ERROR: G4VisCommandReviewKeptEvents::SetNewValue: No kept events,"
	"\n  or kept events not accessible."
	     << G4endl;
    }
    return;
  }

  G4VViewer* viewer = fpVisManager->GetCurrentViewer();
  if (!viewer) {
    if (verbosity >= G4VisManager::errors) {
      G4cout <<
  "ERROR: No current viewer - \"/vis/viewer/list\" to see possibilities."
             << G4endl;
    }
    return;
  }

  G4Scene* pScene = fpVisManager->GetCurrentScene();
  if (!pScene) {
    if (verbosity >= G4VisManager::errors) {
      G4cout << "ERROR: No current scene.  Please create one." << G4endl;
    }
    return;
  }

  G4VSceneHandler* sceneHandler = fpVisManager->GetCurrentSceneHandler();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4int keepVerbose = UImanager->GetVerboseLevel();
  G4int newVerbose(0);
  if (keepVerbose >= 2 || verbosity >= G4VisManager::confirmations)
    newVerbose = 2;
  UImanager->SetVerboseLevel(newVerbose);

  // Event by event refreshing...
  G4bool currentRefreshAtEndOfEvent = pScene->GetRefreshAtEndOfEvent();
  pScene->SetRefreshAtEndOfEvent(true);
  if (macroFileName.empty()) {

    // Draw to viewer and pause session...
    for (size_t i = 0; i < nKeptEvents; ++i) {
      const G4Event* event = (*events)[i];
      if (verbosity >= G4VisManager::warnings) {
	G4cout << "Drawing event : " << event->GetEventID() <<
	  ".  At EndOfEvent, enter any command, then \"continue\"..."
	       << G4endl;
      }
      sceneHandler->SetEvent(event);
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
      if (!viewer->GetViewParameters().IsAutoRefresh())
	UImanager->ApplyCommand("/vis/viewer/flush");
      G4UIsession* session = UImanager->GetSession();
      session->PauseSessionStart("EndOfEvent");
      sceneHandler->SetEvent(0);
    }

  } else {

    // Execute macro file...
    for (size_t i = 0; i < nKeptEvents; ++i) {
      const G4Event* event = (*events)[i];
      if (verbosity >= G4VisManager::warnings) {
	G4cout << "Drawing event : " << event->GetEventID()
	       << " with macro file \"" << macroFileName << G4endl;
      }
      sceneHandler->SetEvent(event);
      UImanager->ApplyCommand("/control/execute " + macroFileName);
      sceneHandler->SetEvent(0);
    }
  }
  pScene->SetRefreshAtEndOfEvent(currentRefreshAtEndOfEvent);

  UImanager->SetVerboseLevel(keepVerbose);
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
