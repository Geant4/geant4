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
// $Id: G4VisCommandsRefresh.hh,v 1.4.2.1 2001/06/28 19:16:09 gunter Exp $
// GEANT4 tag $Name:  $
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
    G4VisManager::PrintCommandDeprecation("Use \"/vis/viewer/refresh\".");
    G4VisManager* pVMan = G4VisManager::GetInstance ();
    if (pVMan -> IsValidView ()) {
      pVMan -> RefreshCurrentView ();
      // Soft clear - clears back buffer only on double-buffered systems.
    }
  }
};

#endif
