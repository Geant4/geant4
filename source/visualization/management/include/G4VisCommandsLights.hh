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
// $Id: G4VisCommandsLights.hh,v 1.6.2.1 2001/06/28 19:16:09 gunter Exp $
// GEANT4 tag $Name:  $
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
    G4VisManager::PrintCommandDeprecation
      ("Use \"/vis/viewer/set/lightsMove\".");
    G4VisManager* pVMan = G4VisManager::GetInstance ();
    pVMan -> SetCurrentViewParameters ().SetLightsMoveWithCamera (value);
    G4VViewer* pView = pVMan -> GetCurrentViewer ();
    if (pView) {
      // Copy current view parameters into current view.
      pView -> SetViewParameters (pVMan -> GetCurrentViewParameters ());
    }
    G4cout << "Issue Draw or refresh to see effect." << G4endl;
  }
};

#endif
