// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsDelete.hh,v 1.4 1999-12-15 14:54:21 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// /vis~/delete/ commands
// John Allison  7th September 1997

#ifndef G4VISCOMMANDSDELETE_HH
#define G4VISCOMMANDSDELETE_HH

#include "globals.hh"
#include "G4VisManager.hh"
#include "G4ViewParameters.hh"

////////////////////////////////////////////////////  /vis~/delete/...  ////
//vis \hline
//vis /vis~/delete/ &&
//vis ...menu of deleting possibilities. \\%
class G4VisCommandDelete {
public:
  G4String GetCommandName () const {return "/vis~/delete/";}
  G4String GetGuidance () const {
    return "...menu of deleting possibilities.";
  }
};

///////////////////////////////////////////  /vis~/delete/scene  ////
//delete \hline
//delete /vis~/delete/scene &&
//delete Deletes the current scene and its views. \\%
class G4VisCommandDeleteScene {
public:
  G4String GetCommandName () const {return "/vis~/delete/scene";}
  G4String GetGuidance () const {
    return "Deletes the current scene and its views.";
  }
  void SetValue () {
    G4VisManager* pVMan = G4VisManager::GetInstance ();
    if (pVMan -> GetVerboseLevel () > 1) {
      pVMan -> PrintCurrentView ();
    }
    pVMan -> DeleteCurrentSceneHandler ();
    if (pVMan -> GetVerboseLevel () > 1) {
      pVMan -> PrintCurrentView ();
    }
  }
};

///////////////////////////////////////////  /vis~/delete/view  ////
//delete \hline
//delete /vis~/delete/view &&
//delete Deletes the current view. \\%
class G4VisCommandDeleteView {
public:
  G4String GetCommandName () const {return "/vis~/delete/view";}
  G4String GetGuidance () const {
    return "Deletes the current view.";
  }
  void SetValue () {
    G4VisManager* pVMan = G4VisManager::GetInstance ();
    if (pVMan -> GetVerboseLevel () > 1) {
      pVMan -> PrintCurrentView ();
    }
    pVMan -> DeleteCurrentViewer ();
    if (pVMan -> GetVerboseLevel () > 1) {
      pVMan -> PrintCurrentView ();
    }
  }
};

#endif
