// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsCreateView.hh,v 1.5 2001-02-05 02:33:53 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// /vis~/create_view/ commands
// John Allison  7th September 1997

#ifndef G4VISCOMMANDSCREATEVIEW_HH
#define G4VISCOMMANDSCREATEVIEW_HH

#include "globals.hh"
#include "G4VisManager.hh"

/////////////////////////////////////////////  /vis~/create_view/...  ////
//vis \hline
//vis /vis~/create\_view/ &&
//vis ...menu of view creation commands. \\%
class G4VisCommandCreateView {
public:
  G4String GetCommandName () const {return "/vis~/create_view/";}
  G4String GetGuidance () const {
    return "...menu of view creation commands.";
  }
};

//  ////////////////////////////////  /vis~/create_view/new_scene ////
//  //cr_vw \hline
//  //cr_vw /vis~/create\_view/new\_scene &&
//  //cr_vw Creates a new scene and a new view; both become current. \\%
class G4VisCommandCreateViewNewScene {
public:
  G4String GetCommandName () const {return "/vis~/create_view/new_scene";}
  G4String GetGuidance () const {
    return "Creates a new scene and a new view; both become current.";
  }
  void SetValue () {
    G4VisManager::PrintCommandDeprecation("Use \"/vis/viewer/create\".");
    G4VisManager* pVMan = G4VisManager::GetInstance ();
    if (pVMan -> IsValidView ()) {
      pVMan -> CreateSceneHandler ();
      pVMan -> CreateViewer ();
      if (pVMan -> GetVerboseLevel () > 1) {
	pVMan -> PrintCurrentScene ();
      }
    }
  }
};

////////////////////////////////  /vis~/create_view/new_view ////
//cr_vw \hline
//cr_vw /vis~/create\_view/new\_view &&
//cr_vw Creates a new view of current scene; new view becomes
//cr_vw current view. \\%
class G4VisCommandCreateViewNewView {
public:
  G4String GetCommandName () const {return "/vis~/create_view/new_view";}
  G4String GetGuidance () const {
    return "Creates a new view of current scene; new view becomes "
      "current view.";
  }
  void SetValue () {
    G4VisManager* pVMan = G4VisManager::GetInstance ();
    if (pVMan -> IsValidView ()) {
      pVMan -> CreateViewer ();
      if (pVMan -> GetVerboseLevel () > 1) {
	pVMan -> PrintCurrentView ();
      }
    }
  }
};

#endif
