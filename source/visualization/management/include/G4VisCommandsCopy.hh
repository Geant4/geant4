// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsCopy.hh,v 1.5 2001-02-05 02:33:51 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// /vis~/copy/ commands
// John Allison  7th September 1997

#ifndef G4VISCOMMANDSCOPY_HH
#define G4VISCOMMANDSCOPY_HH

#include "globals.hh"
#include "G4VisManager.hh"
#include "G4ViewParameters.hh"

////////////////////////////////////////////////////  /vis~/copy/...  ////
//vis \hline
//vis /vis~/copy/ &&
//vis ...menu of copy commands. \\%
class G4VisCommandCopy {
public:
  G4String GetCommandName () const {return "/vis~/copy/";}
  G4String GetGuidance () const {
    return "...menu of copy commands.";
  }
};

/*****************************
//////////////////////////////////////////////  /vis~/copy/all  ////
//copy  \hline
//copy /vis~/copy/all &&
//copy Copy all data of current scene and view into current memory.
//copy WARNING: this overwrites the current scene data and
//copy current view parameters. \\%
class G4VisCommandCopyAll {
public:
  G4String GetCommandName () const {return "/vis~/copy/all";}
  G4String GetGuidance () const {
    return "Copy all data of current scene and view into current memory."
      "  WARNING: this overwrites the current scene data and "
      "current view parameters.";
  }
  void SetValue () {
    G4VisManager* pVMan = G4VisManager::GetInstance ();
    if (pVMan -> IsValidView ()) {
      pVMan -> CopyScene ();
      pVMan -> CopyViewParameters ();
      if (pVMan -> GetVerboseLevel () > 1) {
	pVMan -> PrintCurrentView ();
      }
      if (pVMan -> GetVerboseLevel () > 1) {
	pVMan -> PrintCurrentScene ();
      }
    }
  }
};

//////////////////////////////////////////  /vis~/copy/scene  ////
//copy  \hline
//copy /vis~/copy/scene &&
//copy Copy scene data of current scene into current scene data.
//copy WARNING: this overwrites the current scene data. \\%
class G4VisCommandCopyScene {
public:
  G4String GetCommandName () const {return "/vis~/copy/scene";}
  G4String GetGuidance () const {
    return "Copy scene data of current scene into current scene data."
      "  WARNING: this overwrites the current scene data.";
  }
  void SetValue () {
    G4VisManager* pVMan = G4VisManager::GetInstance ();
    if (pVMan -> IsValidView ()) {
      pVMan -> CopyScene ();
      if (pVMan -> GetVerboseLevel () > 1) {
	pVMan -> PrintCurrentScene ();
      }
    }
  }
};
*********************/

//////////////////////////////////////////  /vis~/copy/view  ////
//copy  \hline
//copy /vis~/copy/view &&
//copy Copy view parameters of current view into current view parameters.
//copy WARNING: this overwrites the current view parameters. \\%
class G4VisCommandCopyView {
public:
  G4String GetCommandName () const {return "/vis~/copy/view";}
  G4String GetGuidance () const {
    return "Copy view parameters of current view into current view parameters."
      "  WARNING: this overwrites the current view parameters.";
  }
  void SetValue () {
    G4VisManager::PrintCommandDeprecation("Use \"/vis/viewer/set/all\".");
    G4VisManager* pVMan = G4VisManager::GetInstance ();
    if (pVMan -> IsValidView ()) {
      pVMan -> CopyViewParameters ();
      if (pVMan -> GetVerboseLevel () > 1) {
	pVMan -> PrintCurrentView ();
      }
    }
  }
};

#endif
