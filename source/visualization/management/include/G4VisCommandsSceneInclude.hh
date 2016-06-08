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
// $Id: G4VisCommandsSceneInclude.hh,v 1.3.4.1 2001/06/28 19:16:10 gunter Exp $
// GEANT4 tag $Name:  $

// /vis/scene commands - John Allison  9th August 1998

#ifndef G4VISCOMMANDSSCENEINCLUDE_HH
#define G4VISCOMMANDSSCENEINCLUDE_HH

#include "G4VisCommandsScene.hh"

class G4VisCommandSceneIncludeHits: public G4VVisCommandScene {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandSceneIncludeHits ();
  ~G4VisCommandSceneIncludeHits ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcmdWithoutParameter* fpCommand;
};

class G4VisCommandSceneIncludeTrajectories: public G4VVisCommandScene {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandSceneIncludeTrajectories ();
  ~G4VisCommandSceneIncludeTrajectories ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcmdWithoutParameter* fpCommand;
};

#endif
