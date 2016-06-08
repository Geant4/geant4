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
// $Id: G4VisCommandsSceneInclude.cc,v 1.7.2.1 2001/06/28 19:16:15 gunter Exp $
// GEANT4 tag $Name:  $

// /vis/scene commands - John Allison  9th August 1998

#include "G4VisCommandsSceneInclude.hh"

#include "G4VisManager.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4HitsModel.hh"
#include "G4TrajectoriesModel.hh"

////////////// /vis/scene/include/hits ///////////////////////////////////////

G4VisCommandSceneIncludeHits::G4VisCommandSceneIncludeHits () {
  fpCommand = new G4UIcmdWithoutParameter ("/vis/scene/include/hits", this);
  fpCommand -> SetGuidance
    ("Replaced by /vis/scene/add/hits.");
}

G4VisCommandSceneIncludeHits::~G4VisCommandSceneIncludeHits () {
  delete fpCommand;
}

G4String G4VisCommandSceneIncludeHits::GetCurrentValue (G4UIcommand* command) {
  return "";
}

void G4VisCommandSceneIncludeHits::SetNewValue (G4UIcommand* command,
						G4String newValue) {
  G4cout <<
    " **** /vis/scene/include/hits has been replaced - \n"
    "  please use the equivalent command /vis/scene/add/hits."
	 << G4endl;
}

////////////// /vis/scene/include/trajectories ///////////////////////////////////////

G4VisCommandSceneIncludeTrajectories::G4VisCommandSceneIncludeTrajectories () {
  fpCommand = new G4UIcmdWithoutParameter
    ("/vis/scene/include/trajectories", this);
  fpCommand -> SetGuidance
    ("Replaced by /vis/scene/add/trajectories.");
}

G4VisCommandSceneIncludeTrajectories::~G4VisCommandSceneIncludeTrajectories () {
  delete fpCommand;
}

G4String G4VisCommandSceneIncludeTrajectories::GetCurrentValue (G4UIcommand* command) {
  return "";
}

void G4VisCommandSceneIncludeTrajectories::SetNewValue (G4UIcommand* command,
					      G4String newValue) {
  G4cout <<
    " **** /vis/scene/include/trajectories has been replaced - \n"
    "  please use the equivalent command /vis/scene/add/trajectories."
	 << G4endl;
}
