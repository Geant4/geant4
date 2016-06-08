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
// $Id: G4VisCommandsCreateScene.hh,v 1.3.4.1 2001/06/28 19:16:09 gunter Exp $
// GEANT4 tag $Name:  $
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
