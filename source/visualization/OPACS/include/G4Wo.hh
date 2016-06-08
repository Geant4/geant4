// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Wo.hh,v 1.5 1999/12/15 14:54:02 gunter Exp $
// GEANT4 tag $Name: geant4-03-00 $
//
// 
// Guy Barrand 04 November 1996
// Wo graphics system factory.

#ifndef G4WO_HH
#define G4WO_HH

#if defined(G4VIS_BUILD_OPACS_DRIVER) || defined(G4VIS_USE_OPACS)

//G4
#include "G4VGraphicsSystem.hh"

class G4VInteractorManager;

class G4Wo: public G4VGraphicsSystem {
public:
  G4Wo ();
  virtual ~G4Wo ();
  G4VSceneHandler*    CreateSceneHandler (const G4String& name = "");
  G4VViewer*     CreateViewer  (G4VSceneHandler&, const G4String& name = "");
  static G4VInteractorManager* GetInteractorManager ();
};

#endif

#endif
