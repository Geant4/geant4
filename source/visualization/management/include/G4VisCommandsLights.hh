// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsLights.hh,v 1.2 1999-01-09 16:30:59 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// /vis~/lights/ commands
// John Allison  7th September 1997

#ifndef G4VISCOMMANDSLIGHTS_HH
#define G4VISCOMMANDSLIGHTS_HH

#include "globals.hh"
#include "G4VisManager.hh"
#include "G4ViewParameters.hh"

////////////////////////////////////////////////////  /vis~/lights/...  ////
//vis \hline
//vis /vis~/lights/ &&
//vis ...menu of lights commands. \\%
class G4VisCommandLights {
public:
  G4String GetCommandName () const {return "/vis~/lights/";}
  G4String GetGuidance () const {
    return "...menu of lights commands.";
  }
};

///////////////////////////////////////  /vis~/lights/move_with_camera ////
//lights \hline
//lights /vis~/lights/move\_with\_camera & true/false &
//lights Lights move with change of viewpoint. \\%
class G4VisCommandLightsMoveWithCamera {
public:
  G4String GetCommandName () const {return "/vis~/lights/move_with_camera";}
  G4String GetGuidance () const {
    return "Lights move with change of viewpoint.";
  }
  G4String GetValueName () const {return "lights move flag";}
  G4bool GetValue () const {
    return G4VisManager::GetInstance () -> GetCurrentViewParameters ().
      GetLightsMoveWithCamera ();
  }
  void SetValue (G4bool value) {
    G4VisManager* pVMan = G4VisManager::GetInstance ();
    pVMan -> SetCurrentViewParameters ().SetLightsMoveWithCamera (value);
    G4VViewer* pView = pVMan -> GetCurrentView ();
    if (pView) {
      // Copy current view parameters into current view.
      pView -> SetViewParameters (pVMan -> GetCurrentViewParameters ());
      // Recalculate projection matrices, etc.
      pView -> SetView ();
    }
    G4cout << "Issue Draw or refresh to see effect." << endl;
  }
};

#endif
