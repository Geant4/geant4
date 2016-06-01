// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4WoView.hh,v 2.1 1998/11/06 13:42:02 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// Guy Barrand 04 November 1996
// Wo view - opens window, hard copy, etc.

#ifndef G4WOVIEW_HH
#define G4WOVIEW_HH

#if defined(G4VIS_BUILD_OPACS_DRIVER) || defined(G4VIS_USE_OPACS)

#include <X11/Intrinsic.h>
#include <OCamera.h>

#include "G4VView.hh"
#include "globals.hh"

class G4GoScene;

// Base class for various WoView classes.
class G4WoView: public G4VView {
public:
          G4WoView  (G4GoScene& scene, const G4String& name = "");
         ~G4WoView  ();
  void    DrawView  ();
  void    ShowView  ();
  OCamera GetCamera ();
private:
  void       ClearView  ();
  void       FinishView ();
  void       SetView    ();
private:
  G4GoScene& fScene;    // Graphics Scene for this view.
  OCamera    fGoCamera;
  Widget     fXoCamera;
  Widget     fShell;
};

#endif

#endif
