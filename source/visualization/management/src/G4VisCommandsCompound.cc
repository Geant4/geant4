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
// $Id: G4VisCommandsCompound.cc 98766 2016-08-09 14:17:17Z gcosmo $

// Compound /vis/ commands - John Allison  15th May 2000

#include "G4VisCommandsCompound.hh"

#include "G4VisManager.hh"
#include "G4UImanager.hh"
#include "G4UIcmdWithAString.hh"

#include <sstream>
#include <set>

////////////// /vis/drawTree ///////////////////////////////////////

G4VisCommandDrawTree::G4VisCommandDrawTree() {
  G4bool omitable;
  fpCommand = new G4UIcommand("/vis/drawTree", this);
  fpCommand->SetGuidance
    ("(DTREE) Creates a scene consisting of this physical volume and"
     "\n  produces a representation of the geometry hieracrhy.");
  fpCommand->SetGuidance("The scene becomes current.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter("physical-volume-name", 's', omitable = true);
  parameter -> SetDefaultValue("world");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter("system", 's', omitable = true);
  parameter -> SetDefaultValue("ATree");
  fpCommand -> SetParameter (parameter);
}

G4VisCommandDrawTree::~G4VisCommandDrawTree() {
  delete fpCommand;
}

void G4VisCommandDrawTree::SetNewValue(G4UIcommand*, G4String newValue) {

  G4String pvname, system;
  std::istringstream is(newValue);
  is >> pvname >> system;

  // Note: The second parameter, "system", is intended to allow the user
  // a choice of dedicated tree printing/displaying systems but at present
  // the only such dedicated system is ASCIITree.  It doesn't make sense to
  // specify OGLSX, for example.  So to avoid confusion we restrict this
  // feature to systems that have "Tree" in the name or nickname.

  // Of course, some other systems, such as OGLSQt, have a tree browser
  // built-in.  The HepRApp offline browser also has a tree browser
  // built in.

  if (!system.contains("Tree")) {
    system = "ATree";
  }

  G4VGraphicsSystem* keepSystem = fpVisManager->GetCurrentGraphicsSystem();
  G4Scene* keepScene = fpVisManager->GetCurrentScene();
  G4VSceneHandler* keepSceneHandler = fpVisManager->GetCurrentSceneHandler();
  G4VViewer* keepViewer = fpVisManager->GetCurrentViewer();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4int keepVerbose = UImanager->GetVerboseLevel();
  G4int newVerbose(0);
  if (keepVerbose >= 2 ||
      fpVisManager->GetVerbosity() >= G4VisManager::confirmations)
    newVerbose = 2;
  UImanager->SetVerboseLevel(newVerbose);
  UImanager->ApplyCommand(G4String("/vis/open " + system));
  if (fErrorCode == 0) {
    UImanager->ApplyCommand(G4String("/vis/drawVolume " + pvname));
    UImanager->ApplyCommand("/vis/viewer/flush");
    if (keepViewer) {
      if (fpVisManager->GetVerbosity() >= G4VisManager::warnings) {
        G4cout << "Reverting to " << keepViewer->GetName() << G4endl;
      }
      fpVisManager->SetCurrentGraphicsSystem(keepSystem);
      fpVisManager->SetCurrentScene(keepScene);
      fpVisManager->SetCurrentSceneHandler(keepSceneHandler);
      fpVisManager->SetCurrentViewer(keepViewer);
    }
  }
  UImanager->SetVerboseLevel(keepVerbose);
}

////////////// /vis/drawView ///////////////////////////////////////

G4VisCommandDrawView::G4VisCommandDrawView() {
  G4bool omitable;
  fpCommand = new G4UIcommand("/vis/drawView", this);
  fpCommand->SetGuidance
    ("Draw view from this angle, etc.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter("theta-degrees", 'd', omitable = true);
  parameter -> SetDefaultValue(0.);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter("phi-degrees", 'd', omitable = true);
  parameter -> SetDefaultValue(0.);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter("pan-right", 'd', omitable = true);
  parameter -> SetDefaultValue(0.);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter("pan-up", 'd', omitable = true);
  parameter -> SetDefaultValue(0.);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter("pan-unit", 's', omitable = true);
  parameter -> SetDefaultValue("cm");
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter("zoom-factor", 'd', omitable = true);
  parameter -> SetDefaultValue(1.);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter("dolly", 'd', omitable = true);
  parameter -> SetDefaultValue(0.);
  fpCommand -> SetParameter (parameter);
  parameter = new G4UIparameter("dolly-unit", 's', omitable = true);
  parameter -> SetDefaultValue("cm");
  fpCommand -> SetParameter (parameter);
}

G4VisCommandDrawView::~G4VisCommandDrawView() {
  delete fpCommand;
}

void G4VisCommandDrawView::SetNewValue(G4UIcommand*, G4String newValue) {

  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();

  G4VViewer* currentViewer = fpVisManager->GetCurrentViewer();
  if (!currentViewer) {
    if (verbosity >= G4VisManager::warnings) {
      G4cout <<
	"WARNING: G4VisCommandsDrawView::SetNewValue: no current viewer."
	     << G4endl;
    }
    return;
  }

  G4String thetaDeg;
  G4String phiDeg;
  G4String panRight;
  G4String panUp;
  G4String panUnit;
  G4String zoomFactor;
  G4String dolly;
  G4String dollyUnit;
  std::istringstream is(newValue);
  is >> thetaDeg >> phiDeg >> panRight >> panUp >> panUnit
     >> zoomFactor >> dolly >> dollyUnit;
  
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4int keepVerbose = UImanager->GetVerboseLevel();
  G4int newVerbose(0);
  if (keepVerbose >= 2 ||
      fpVisManager->GetVerbosity() >= G4VisManager::confirmations)
    newVerbose = 2;
  UImanager->SetVerboseLevel(newVerbose);
  G4ViewParameters vp = currentViewer->GetViewParameters();
  G4bool keepAutoRefresh = vp.IsAutoRefresh();
  vp.SetAutoRefresh(false);
  currentViewer->SetViewParameters(vp);
  UImanager->ApplyCommand(
    G4String("/vis/viewer/set/viewpointThetaPhi " + thetaDeg + " " + phiDeg + " deg"));
  UImanager->ApplyCommand(
    G4String("/vis/viewer/panTo " + panRight + " " + panUp + " " + panUnit));
  UImanager->ApplyCommand(
    G4String("/vis/viewer/zoomTo " + zoomFactor));
  vp = currentViewer->GetViewParameters();
  vp.SetAutoRefresh(keepAutoRefresh);
  currentViewer->SetViewParameters(vp);
  UImanager->ApplyCommand(
    G4String("/vis/viewer/dollyTo " + dolly + " " + dollyUnit));
  UImanager->SetVerboseLevel(keepVerbose);
}

////////////// /vis/drawVolume ///////////////////////////////////////

G4VisCommandDrawVolume::G4VisCommandDrawVolume() {
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString("/vis/drawVolume", this);
  fpCommand->SetGuidance
    ("Creates a scene containing this physical volume and asks the"
     "\ncurrent viewer to draw it.  The scene becomes current.");
  fpCommand -> SetGuidance 
    ("If physical-volume-name is \"world\" (the default), the main geometry"
     "\ntree (material world) is drawn.  If \"worlds\", all worlds - material"
     "\nworld and parallel worlds, if any - are drawn.  Otherwise a search of"
     "\nall worlds is made, taking the first matching occurence only.  To see"
     "\na representation of the geometry hierarchy of the worlds, try"
     "\n\"/vis/drawTree worlds\" or one of the driver/browser combinations"
     "\nthat have the required functionality, e.g., HepRepFile/HepRApp.");
  fpCommand->SetParameterName("physical-volume-name", omitable = true);
  fpCommand->SetDefaultValue("world");
}

G4VisCommandDrawVolume::~G4VisCommandDrawVolume() {
  delete fpCommand;
}

void G4VisCommandDrawVolume::SetNewValue(G4UIcommand*, G4String newValue) {
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4int keepVerbose = UImanager->GetVerboseLevel();
  G4int newVerbose(0);
  if (keepVerbose >= 2 || verbosity >= G4VisManager::confirmations)
    newVerbose = 2;
  UImanager->SetVerboseLevel(newVerbose);
  UImanager->ApplyCommand("/vis/scene/create");
  UImanager->ApplyCommand(G4String("/vis/scene/add/volume " + newValue));
  UImanager->ApplyCommand("/vis/sceneHandler/attach");
  UImanager->SetVerboseLevel(keepVerbose);
  static G4bool warned = false;
  if (verbosity >= G4VisManager::confirmations && !warned) {
    G4cout <<
      "NOTE: For systems which are not \"auto-refresh\" you will need to"
      "\n  issue \"/vis/viewer/refresh\" or \"/vis/viewer/flush\"."
	   << G4endl;
    warned = true;
  }
}

////////////// /vis/open ///////////////////////////////////////

G4VisCommandOpen::G4VisCommandOpen() {
  G4bool omitable;
  fpCommand = new G4UIcommand("/vis/open", this);
  fpCommand->SetGuidance
    ("Creates a scene handler ready for drawing.");
  fpCommand->SetGuidance
    ("The scene handler becomes current (the name is auto-generated).");
  G4UIparameter* parameter;
  parameter = new G4UIparameter("graphics-system-name", 's', omitable = false);
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("window-size-hint", 's', omitable = true);
  parameter->SetGuidance
    ("integer (pixels) for square window placed by window manager or"
     " X-Windows-type geometry string, e.g. 600x600-100+100");
  parameter->SetDefaultValue("600");
  fpCommand->SetParameter(parameter);
}

G4VisCommandOpen::~G4VisCommandOpen() {
  delete fpCommand;
}

void G4VisCommandOpen::SetNewValue (G4UIcommand*, G4String newValue) {
  G4String systemName, windowSizeHint;
  std::istringstream is(newValue);
  is >> systemName >> windowSizeHint;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4int keepVerbose = UImanager->GetVerboseLevel();
  G4int newVerbose(0);
  if (keepVerbose >= 2 ||
      fpVisManager->GetVerbosity() >= G4VisManager::confirmations)
    newVerbose = 2;
  UImanager->SetVerboseLevel(newVerbose);
  fErrorCode = UImanager->ApplyCommand(G4String("/vis/sceneHandler/create " + systemName));
  if (fErrorCode == 0) {
    UImanager->ApplyCommand(G4String("/vis/viewer/create ! ! " + windowSizeHint));
  } else {
    // Use set to get alphabetical order
    std::set<G4String> candidates;
    for (const auto gs: fpVisManager -> GetAvailableGraphicsSystems()) {
      // Just list nicknames, but exclude FALLBACK nicknames
      for (const auto& nickname: gs->GetNicknames()) {
        if (!nickname.contains("FALLBACK")) {
          candidates.insert(nickname);
        }
      }
    }
    G4cerr << "Candidates are:";
    for (const auto& candidate: candidates) {
      G4cerr << ' ' << candidate;
    }
    G4cerr << G4endl;
  }
  UImanager->SetVerboseLevel(keepVerbose);
}

////////////// /vis/specify ///////////////////////////////////////

G4VisCommandSpecify::G4VisCommandSpecify() {
  G4bool omitable;
  fpCommand = new G4UIcommand("/vis/specify", this);
  fpCommand->SetGuidance
    ("Draws logical volume with Boolean components, voxels and readout geometry.");
  fpCommand->SetGuidance
    ("Creates a scene consisting of this logical volume and asks the"
     "\n  current viewer to draw it to the specified depth of descent"
     "\n  showing boolean components (if any), voxels (if any),"
     "\n  readout geometry (if any) and local axes, under control of the"
     "\n  appropriate flag.");
  fpCommand->SetGuidance
  ("Note: voxels are not constructed until start of run - /run/beamOn."
   "\n  (For voxels without a run, \"/run/beamOn 0\".)");
  fpCommand->SetGuidance("The scene becomes current.");
  G4UIparameter* parameter;
  parameter = new G4UIparameter("logical-volume-name", 's', omitable = false);
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("depth-of-descent", 'i', omitable = true);
  parameter->SetDefaultValue(1);
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("booleans-flag", 'b', omitable = true);
  parameter->SetDefaultValue(true);
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("voxels-flag", 'b', omitable = true);
  parameter->SetDefaultValue(true);
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("readout-flag", 'b', omitable = true);
  parameter->SetDefaultValue(true);
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("axes-flag", 'b', omitable = true);
  parameter->SetDefaultValue(true);
  parameter -> SetGuidance ("Set \"false\" to suppress axes.");
  fpCommand->SetParameter(parameter);
}

G4VisCommandSpecify::~G4VisCommandSpecify() {
  delete fpCommand;
}

void G4VisCommandSpecify::SetNewValue(G4UIcommand*, G4String newValue) {
  G4VisManager::Verbosity verbosity = fpVisManager->GetVerbosity();
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4int keepVerbose = UImanager->GetVerboseLevel();
  G4int newVerbose(0);
  if (keepVerbose >= 2 || verbosity >= G4VisManager::confirmations)
    newVerbose = 2;
  UImanager->SetVerboseLevel(newVerbose);
  // UImanager->ApplyCommand(G4String("/geometry/print " + newValue));
  UImanager->ApplyCommand("/vis/scene/create");
  UImanager->ApplyCommand(G4String("/vis/scene/add/logicalVolume " + newValue));
  UImanager->ApplyCommand("/vis/sceneHandler/attach");
  UImanager->SetVerboseLevel(keepVerbose);
  static G4bool warned = false;
  if (verbosity >= G4VisManager::confirmations && !warned) {
    G4cout <<
      "NOTE: For systems which are not \"auto-refresh\" you will need to"
      "\n  issue \"/vis/viewer/refresh\" or \"/vis/viewer/flush\"."
	   << G4endl;
    warned = true;
  }
}
