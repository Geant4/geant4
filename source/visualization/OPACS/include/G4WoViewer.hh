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
// $Id: G4WoViewer.hh,v 1.5 2001/07/11 10:08:47 gunter Exp $
// GEANT4 tag $Name: geant4-05-00 $
//
// 
// Guy Barrand 04 November 1996
// Wo viewer - opens window, hard copy, etc.

#ifndef G4WOVIEWER_HH
#define G4WOVIEWER_HH

#if defined(G4VIS_BUILD_OPACS_DRIVER) || defined(G4VIS_USE_OPACS)

#include <X11/Intrinsic.h>
#include <OCamera.h>

#include "G4VViewer.hh"
#include "globals.hh"

class G4GoSceneHandler;

// Base class for various WoView classes.
class G4WoViewer: public G4VViewer {
public:
          G4WoViewer  (G4GoSceneHandler& scene, const G4String& name = "");
 virtual ~G4WoViewer  ();
  void    DrawView  ();
  void    ShowView  ();
  OCamera GetCamera ();
private:
  void       ClearView  ();
  void       FinishView ();
  void       SetView    ();
private:
  G4GoSceneHandler& fSceneHandler;    // Graphics Scene for this view.
  OCamera    fGoCamera;
  Widget     fXoCamera;
  Widget     fShell;
};

#endif

#endif
