// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsCreateScene.hh,v 1.2 1999-01-09 16:30:57 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// /vis~/create_scene/ commands
// John Allison  7th September 1997

#ifndef G4VISCOMMANDSCREATESCENE_HH
#define G4VISCOMMANDSCREATESCENE_HH

#include "globals.hh"
#include "G4UImanager.hh"

/////////////////////////////////////////////  /vis~/create_scene/...  ////
//vis \hline
//vis /vis~/create\_scene/ &&
//vis ...menu of scene creation commands. \\%
class G4VisCommandCreateScene {
public:
  G4String GetCommandName () const {return "/vis~/create_scene/";}
  G4String GetGuidance () const {
    return "...menu of scene creation commands.";
  }
};

//////////////////////////////////  /vis~/create_scene/clear  ////
//cr_sc \hline
//cr_sc /vis~/create\_scene/clear &&
//cr_sc Same as {\tt /vis~/clear/scene}. \\%
class G4VisCommandCreateSceneClear {
public:
  G4String GetCommandName () const {return "/vis~/create_scene/clear";}
  G4String GetGuidance () const {
    return "Same as /vis~/clear/scene.";
  }
  void SetValue () {
    G4UImanager* UI = G4UImanager::GetUIpointer ();
    UI -> ApplyCommand("/vis~/clear/scene");
  }
};

#endif
