// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsShow.hh,v 1.1 1999-01-07 16:15:23 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// /vis~/show/ commands
// John Allison  7th September 1997

#ifndef G4VISCOMMANDSSHOW_HH
#define G4VISCOMMANDSSHOW_HH

#include "globals.hh"
#include "G4VisManager.hh"

///////////////////////////////////////////////////  /vis~/show/...  ////
//vis \hline
//vis /vis~/show/ &&
//vis ...menu of show commands. \\%
class G4VisCommandShow {
public:
  G4String GetCommandName () const {return "/vis~/show/";}
  G4String GetGuidance () const {
    return "...menu of show commands.";
  }
};

////////////////////////////////////////////  /vis~/show/view  ////
//view \hline
//view /vis~/show/view &&
//view Make the view visible, if not so already (initiates
//view post-processing for
//view graphics systems which use such techniques). \\%
class G4VisCommandShowView {
public:
  G4String GetCommandName () const {return "/vis~/show/view";}
  G4String GetGuidance () const {
    return "Make the view visible, if not so already (initiates "
      "post-processing for "
      "graphics systems which use such techniques).";
  }
  void SetValue () {
    G4VisManager* pVMan = G4VisManager::GetInstance ();
    if (pVMan -> IsValidView ()) {
      if (pVMan -> GetVerboseLevel () > 1) {
	pVMan -> PrintCurrentView ();
      }
      pVMan -> Show ();    
    }
  }
};

#endif
