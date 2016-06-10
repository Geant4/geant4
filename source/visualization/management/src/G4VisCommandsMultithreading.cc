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
// $Id: G4VisCommandsMultithreading.cc 88742 2015-03-08 13:12:57Z allison $

// /vis/multithreading commands - John Allison  29th September 2015

#ifdef G4MULTITHREADED

#include "G4VisCommandsMultithreading.hh"

#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4VisManager.hh"

////////////// /vis/multithreading/actionOnEventQueueFull ///////////////////////////////////////

G4VisCommandMultithreadingActionOnEventQueueFull::G4VisCommandMultithreadingActionOnEventQueueFull()
{
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString("/vis/multithreading/actionOnEventQueueFull", this);
  fpCommand->SetGuidance("When event queue for drawing gets full:");
  fpCommand->SetGuidance("wait: event processing waits for vis manager to catch up.");
  fpCommand->SetGuidance("discard: events are discarded for drawing.");
  fpCommand->SetCandidates("wait discard");
  fpCommand->SetParameterName ("wait", omitable = true);
  fpCommand->SetDefaultValue ("wait");
}

G4VisCommandMultithreadingActionOnEventQueueFull::~G4VisCommandMultithreadingActionOnEventQueueFull()
{
  delete fpCommand;
}

G4String G4VisCommandMultithreadingActionOnEventQueueFull::GetCurrentValue(G4UIcommand*)
{
  return "";
}

void G4VisCommandMultithreadingActionOnEventQueueFull::SetNewValue(G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  if (newValue == "wait") {
    fpVisManager->SetWaitOnEventQueueFull(true);
  } else {
    fpVisManager->SetWaitOnEventQueueFull(false);
  }

  if (verbosity >= G4VisManager::confirmations) {
    G4cout <<
    "When event queue for drawing is full,";
    if (fpVisManager->GetWaitOnEventQueueFull()) {
      G4cout << " event processing will wait";
    } else {
      G4cout << " events will be discarded for drawing";
    }
    G4cout << G4endl;
  }
}

////////////// /vis/multithreading/maxEventQueueSize ///////////////////////////////////////

G4VisCommandMultithreadingMaxEventQueueSize::G4VisCommandMultithreadingMaxEventQueueSize()
{
  G4bool omitable;
  fpCommand = new G4UIcmdWithAnInteger("/vis/multithreading/maxEventQueueSize", this);
  fpCommand->SetGuidance
  ("Defines maximum event queue size. N <=0 means \"unlimited\".");
  fpCommand->SetGuidance
  ("If adding an event to the visualisation event queue would cause the"
   " queue size to exceed this value:");
  fpCommand->SetGuidance
  (" if actionOnEventQueueFull==wait the worker threads are paused for a short"
   " time to give the visualisation manager a chance to catch up.");
  fpCommand->SetGuidance
  (" if actionOnEventQueueFull==discard the event is discarded for drawing.");
  fpCommand->SetParameterName ("maxSize", omitable = true);
  fpCommand->SetDefaultValue (100);
}

G4VisCommandMultithreadingMaxEventQueueSize::~G4VisCommandMultithreadingMaxEventQueueSize()
{
  delete fpCommand;
}

G4String G4VisCommandMultithreadingMaxEventQueueSize::GetCurrentValue(G4UIcommand*)
{
  return "";
}

void G4VisCommandMultithreadingMaxEventQueueSize::SetNewValue(G4UIcommand*, G4String newValue)
{
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4int maxEventQueueSize = fpCommand->GetNewIntValue(newValue);
  fpVisManager->SetMaxEventQueueSize(maxEventQueueSize);

  if (verbosity >= G4VisManager::confirmations) {
    G4cout <<
    "Maximum event queue size has been set to "
	   << fpVisManager->GetMaxEventQueueSize()
	   << G4endl;
  }
}

#endif
