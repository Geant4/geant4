// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsCompound.cc,v 1.3 2000-06-02 12:24:53 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// Compound /vis/ commands - John Allison  15th May 2000

#include "G4VisCommandsCompound.hh"

#include "G4VisManager.hh"
#include "G4UImanager.hh"
#include "G4UIcmdWithAString.hh"

////////////// /vis/drawVolume ///////////////////////////////////////

G4VisCommandDrawVolume::G4VisCommandDrawVolume() {
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString("/vis/drawVolume", this);
  fpCommand->AvailableForStates(Idle, GeomClosed);
  fpCommand->SetGuidance("/vis/drawVolume [<physical-volume-name>]");
  fpCommand->SetGuidance("Default: world volume");
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

void G4VisCommandDrawVolume::SetNewValue
(G4UIcommand* command, G4String newValue) {
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4int keepVerbose = UImanager->GetVerboseLevel();
  UImanager->SetVerboseLevel(2);
  UImanager->ApplyCommand("/vis/scene/create");
  UImanager->ApplyCommand("/vis/scene/add/volume " + newValue);
  UImanager->ApplyCommand("/vis/sceneHandler/attach");
  UImanager->ApplyCommand("/vis/viewer/refresh");
  UImanager->ApplyCommand("/vis/viewer/show");
  UImanager->SetVerboseLevel(keepVerbose);
}

////////////// /vis/open ///////////////////////////////////////

G4VisCommandOpen::G4VisCommandOpen() {
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString("/vis/open", this);
  fpCommand->SetGuidance("/vis/open <graphics-system-name>");
  fpCommand->SetGuidance
    ("For this graphics system, creates a scene handler ready for drawing.");
  fpCommand->SetGuidance("The scene handler becomes current.");
  fpCommand->SetGuidance("The scene handler name is auto-generated.");
  fpCommand->SetParameterName("graphics-system-name", omitable = false);
  const G4GraphicsSystemList& gslist =
    fpVisManager->GetAvailableGraphicsSystems();
  G4String candidates;
  for (int igslist = 0; igslist < gslist.entries(); igslist++) {
    const G4String& name = gslist(igslist)->GetName();
    const G4String& nickname = gslist(igslist)->GetNickname();
    if (nickname.isNull()) {
      candidates += name;
    }
    else {
      candidates += nickname;
    }
    candidates += " ";
  }
  candidates = candidates.strip();
  fpCommand->SetCandidates(candidates);
}

G4VisCommandOpen::~G4VisCommandOpen() {
  delete fpCommand;
}

void G4VisCommandOpen::SetNewValue (G4UIcommand* command, G4String newValue) {
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/vis/sceneHandler/create " + newValue);
  UImanager->ApplyCommand("/vis/viewer/create");
}

////////////// /vis/specify ///////////////////////////////////////

G4VisCommandSpecify::G4VisCommandSpecify() {
  G4bool omitable;
  fpCommand = new G4UIcmdWithAString("/vis/specify", this);
  fpCommand->AvailableForStates(Idle, GeomClosed);
  fpCommand->SetGuidance("/vis/specify <logical-volume-name>");
  fpCommand->SetGuidance
    ("Creates a scene consisting of this logical volume and asks the"
     "\n  current viewer to draw it and the geometry to print the"
     "\n  specification.");
  fpCommand->SetGuidance("The scene becomes current.");
  fpCommand->SetParameterName("logical-volume-name", omitable = false);
}

G4VisCommandSpecify::~G4VisCommandSpecify() {
  delete fpCommand;
}

void G4VisCommandSpecify::SetNewValue
(G4UIcommand* command, G4String newValue) {
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4int keepVerbose = UImanager->GetVerboseLevel();
  UImanager->SetVerboseLevel(2);
  UImanager->ApplyCommand("/vis/scene/create");
  UImanager->ApplyCommand("/vis/scene/add/logicalVolume " + newValue);
  UImanager->ApplyCommand("/vis/sceneHandler/attach");
  UImanager->ApplyCommand("/vis/viewer/refresh");
  UImanager->ApplyCommand("/vis/viewer/show");
  UImanager->ApplyCommand("/geometry/print " + newValue);
  UImanager->SetVerboseLevel(keepVerbose);
}
