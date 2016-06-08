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
// $Id: G4VisCommandsClear.hh,v 1.5.2.1 2001/06/28 19:16:08 gunter Exp $
// GEANT4 tag $Name:  $
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
    G4VisManager::PrintCommandDeprecation("Use \"/vis/viewer/clear\".");
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
