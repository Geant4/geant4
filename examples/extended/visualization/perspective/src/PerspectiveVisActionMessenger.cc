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
/// \file visualization/perspective/src/PerspectiveVisActionMessenger.cc
/// \brief Implementation of the PerspectiveVisActionMessenger class
//
//

#include "PerspectiveVisActionMessenger.hh"

#include "PerspectiveVisAction.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PerspectiveVisActionMessenger::PerspectiveVisActionMessenger
(PerspectiveVisAction* PVA):
  G4UImessenger(),
  fPVA(PVA),
  fpDirectory(0),
  fpCommandOS(0),
  fpCommandScene(0)
{
  G4bool omitable;

  fpDirectory = new G4UIdirectory ("/perspectiveDemo/");
  fpDirectory -> SetGuidance ("Perspective demonstration commands.");

  fpCommandOS = new G4UIcmdWithAString ("/perspectiveDemo/optionString", this);
  fpCommandOS -> SetGuidance
    ("Option string - any combination of \"x\", \"y\", \"z\", \"a[ll]\".");
  fpCommandOS -> SetGuidance
    ("Controls direction of perspective lines.");
  fpCommandOS -> SetParameterName ("option-string", omitable = true);
  fpCommandOS -> SetDefaultValue("all");

  fpCommandScene = new G4UIcmdWithAString ("/perspectiveDemo/scene", this);
  fpCommandScene -> SetGuidance
    ("Scene name.");
  fpCommandScene -> SetParameterName ("scene-name", omitable = true);
  fpCommandScene -> SetDefaultValue("room-and-chair");
  fpCommandScene -> SetCandidates("room-and-chair");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PerspectiveVisActionMessenger::~PerspectiveVisActionMessenger ()
{
  delete fpCommandScene;
  delete fpCommandOS;
  delete fpDirectory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

