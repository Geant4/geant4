// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsClear.hh,v 1.3 1999-01-11 00:48:19 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// /vis~/clear/ commands
// John Allison  7th September 1997

#ifndef G4VISCOMMANDSCLEAR_HH
#define G4VISCOMMANDSCLEAR_HH

#include "globals.hh"
#include "G4VisManager.hh"
#include "G4ViewParameters.hh"

////////////////////////////////////////////////////  /vis~/clear/...  ////
//vis \hline
//vis /vis~/clear/ &&
//vis ...menu of clear commands. \\%
class G4VisCommandClear {
public:
  G4String GetCommandName () const {return "/vis~/clear/";}
  G4String GetGuidance () const {
    return "...menu of clear commands.";
  }
};

/****************
/////////////////////////////////////////  /vis~/clear/scene  ////
//clear \hline
//clear /vis~/clear/scene &&
//clear Clears current scene and scene data and marks all its views
//clear as needing refreshing (but does not clear them). \\%
class G4VisCommandClearScene {
public:
  G4String GetCommandName () const {return "/vis~/clear/scene";}
  G4String GetGuidance () const {
    return  "Clears current scene and scene data and marks all its views"
      " as needing refreshing (but does not clear them).";
  }
  void SetValue () {
    G4VisManager* pVMan = G4VisManager::GetInstance ();
    if (pVMan -> IsValidView ()) {
      pVMan -> ClearScene ();    
      if (pVMan -> GetVerboseLevel () > 1) {
	pVMan -> PrintCurrentView ();
      }
    }
  }
};
******************/

////////////////////////////////////////  /vis~/clear/view ////
//clear \hline
//clear /vis~/clear/view &&
//clear Clears visible window of current view. \\%
class G4VisCommandClearView {
public:
  G4String GetCommandName () const {return "/vis~/clear/view";}
  G4String GetGuidance () const {
    return "Clears visible window of current view.";
  }
  void SetValue () {
    G4VisManager* pVMan = G4VisManager::GetInstance ();
    if (pVMan -> IsValidView ()) {
      pVMan -> ClearView ();    
      if (pVMan -> GetVerboseLevel () > 1) {
	pVMan -> PrintCurrentView ();
      }
    }
  }
};

/****************************
///////////////////////////////////////  /vis~/clear/view_and_scene  ////
//clear \hline
//clear /vis~/clear/view\_and\_scene &&
//clear Combination command. \\%
class G4VisCommandClearViewAndScene {
public:
  G4String GetCommandName () const {return "/vis~/clear/view_and_scene";}
  G4String GetGuidance () const {
    return "Combination command.";
  }
  void SetValue () {
    G4VisManager* pVMan = G4VisManager::GetInstance ();
    if (pVMan -> IsValidView ()) {
      pVMan -> ClearView ();    
      if (pVMan -> GetVerboseLevel () > 1) {
	pVMan -> PrintCurrentView ();
      }
      pVMan -> ClearScene ();    
      if (pVMan -> GetVerboseLevel () > 1) {
	pVMan -> PrintCurrentView ();
      }
    }
  }
};
********************/

#endif
