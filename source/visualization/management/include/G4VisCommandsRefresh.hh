// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsRefresh.hh,v 1.2 1999-01-09 16:31:00 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// /vis~/refresh/ commands
// John Allison  7th September 1997

#ifndef G4VISCOMMANDSREFRESH_HH
#define G4VISCOMMANDSREFRESH_HH

#include "globals.hh"
#include "G4VisManager.hh"

///////////////////////////////////////////////////  /vis~/refresh/...  ////
//vis \hline
//vis /vis~/refresh/ &&
//vis ...menu of refresh commands. \\%
class G4VisCommandRefresh {
public:
  G4String GetCommandName () const {return "/vis~/refresh/";}
  G4String GetGuidance () const {
    return "...menu of refresh commands.";
  }
};

/////////////////////////////////////////////  /vis~/refresh/view  ////
//view \hline
//view /vis~/refresh/view &&
//view Optimal redraw of current view of current scene - uses
//view double buffer and graphical database, if any. \\%
class G4VisCommandRefreshView {
public:
  G4String GetCommandName () const {return "/vis~/refresh/view";}
  G4String GetGuidance () const {
    return "Optimal redraw of current view of current scene - uses "
      "double buffer and graphical database, if any.";
  }
  void SetValue () {
    G4VisManager* pVMan = G4VisManager::GetInstance ();
    if (pVMan -> IsValidView ()) {
      pVMan -> RefreshCurrentView ();
      // Soft clear - clears back buffer only on double-buffered systems.
    }
  }
};

#endif
