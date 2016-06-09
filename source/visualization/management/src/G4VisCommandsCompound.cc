//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VisCommandsCompound.cc,v 1.28 2005/03/03 16:13:08 allison Exp $
// GEANT4 tag $Name: geant4-07-01 $

// Compound /vis/ commands - John Allison  15th May 2000

#include "G4VisCommandsCompound.hh"

#include "G4VisManager.hh"
#include "G4UImanager.hh"
#include "G4UIcmdWithAString.hh"

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
  std::istrstream is(newValue);
  is >> pvname >> system;

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4int keepVerbose = UImanager->GetVerboseLevel();
  G4int newVerbose(0);
  if (keepVerbose >= 2 ||
      fpVisManager->GetVerbosity() >= G4VisManager::confirmations)
    newVerbose = 2;
  UImanager->SetVerboseLevel(newVerbose);
  UImanager->ApplyCommand(G4String("/vis/open " + system));
  UImanager->ApplyCommand(G4String("/vis/drawVolume " + pvname));
  UImanager->ApplyCommand("/vis/viewer/flush");
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
  std::istrstream is(newValue);
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
    ("Creates a scene consisting of this physical volume and asks the"
     "\n  current viewer to draw it.");
  fpCommand->SetGuidance("The scene becomes current.");
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
  if (verbosity >= G4VisManager::warnings) {
    G4cout <<
      "WARNING: For systems which are not \"auto-refresh\" you will need to"
      "\n  issue \"/vis/viewer/refresh\" or \"/vis/viewer/flush\"."
	   << G4endl;
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
   const G4GraphicsSystemList& gslist =
    fpVisManager->GetAvailableGraphicsSystems();
  G4String candidates;
  for (size_t igslist = 0; igslist < gslist.size(); igslist++) {
    const G4String& name = gslist[igslist]->GetName();
    const G4String& nickname = gslist[igslist]->GetNickname();
    if (nickname.isNull()) {
      candidates += name;
    }
    else {
      candidates += nickname;
    }
    candidates += " ";
  }
  candidates = candidates.strip();
  parameter->SetParameterCandidates(candidates);
  fpCommand->SetParameter(parameter);
  parameter = new G4UIparameter("window-size-hint", 'i', omitable = true);
  parameter->SetGuidance("pixels");
  parameter->SetDefaultValue(600);
  fpCommand->SetParameter(parameter);
}

G4VisCommandOpen::~G4VisCommandOpen() {
  delete fpCommand;
}

void G4VisCommandOpen::SetNewValue (G4UIcommand*, G4String newValue) {
  G4String systemName, windowSizeHint;
  std::istrstream is(newValue);
  is >> systemName >> windowSizeHint;
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4int keepVerbose = UImanager->GetVerboseLevel();
  G4int newVerbose(0);
  if (keepVerbose >= 2 ||
      fpVisManager->GetVerbosity() >= G4VisManager::confirmations)
    newVerbose = 2;
  UImanager->SetVerboseLevel(newVerbose);
  UImanager->ApplyCommand(G4String("/vis/sceneHandler/create " + systemName));
  UImanager->ApplyCommand(G4String("/vis/viewer/create ! ! " + windowSizeHint));
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
     "\n  showing boolean components (if any), voxels (if any)"
     "\n  and readout geometry (if any), under control of the appropriate flag.");
  fpCommand->SetGuidance
    ("Note: voxels are not constructed until start of run - /run/beamOn.");
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
  if (verbosity >= G4VisManager::warnings) {
    G4cout <<
      "WARNING: For systems which are not \"auto-refresh\" you will need to"
      "\n  issue \"/vis/viewer/refresh\" or \"/vis/viewer/flush\"."
	   << G4endl;
  }
}
