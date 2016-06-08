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
// $Id: G4VisCommandsDelete.hh,v 1.5.2.1 2001/06/28 19:16:09 gunter Exp $
// GEANT4 tag $Name:  $
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
    G4VisManager::PrintCommandDeprecation
      ("Use \"/vis/scene/remove\" and/or \"/vis/sceneHandler/remove\".");
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
    G4VisManager::PrintCommandDeprecation("Use \"/vis/viewer/remove\".");
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
