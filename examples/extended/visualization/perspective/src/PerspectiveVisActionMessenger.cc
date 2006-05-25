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
// $Id: PerspectiveVisActionMessenger.cc,v 1.1 2006-05-25 08:44:52 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "PerspectiveVisActionMessenger.hh"

#include "PerspectiveVisAction.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

PerspectiveVisActionMessenger::PerspectiveVisActionMessenger
(PerspectiveVisAction* PVA):
  fPVA(PVA)
{
  G4bool omitable;

  fpDirectory = new G4UIdirectory ("/perspectiveDemo/");
  fpDirectory -> SetGuidance ("Perspective demonstration commands.");

  fpCommandOS = new G4UIcmdWithAString ("/perspectiveDemo/optionString", this);
  fpCommandOS -> SetGuidance
    ("Option string - passed in G4AttValue."
     "\n Any combination of \"x\", \"y\", \"z\", \"a[ll]\".");
  fpCommandOS -> SetParameterName ("option-string", omitable = true);
  fpCommandOS -> SetDefaultValue("all");

  fpCommandScene = new G4UIcmdWithAString ("/perspectiveDemo/scene", this);
  fpCommandScene -> SetGuidance
    ("Scene name.");
  fpCommandScene -> SetParameterName ("scene-name", omitable = true);
  fpCommandScene -> SetDefaultValue("room-and-chair");
  fpCommandScene -> SetCandidates("room-and-chair");
}

PerspectiveVisActionMessenger::~PerspectiveVisActionMessenger ()
{
  delete fpCommandScene;
  delete fpCommandOS;
  delete fpDirectory;
}

void PerspectiveVisActionMessenger::SetNewValue
(G4UIcommand* command, G4String newValue)
{
  if (command == fpCommandOS)
    {
      fPVA->SetOptionString(newValue);
    }

  else if (command == fpCommandScene)
    {
      fPVA->SetScene(newValue);
    }

  G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/rebuild");
}
