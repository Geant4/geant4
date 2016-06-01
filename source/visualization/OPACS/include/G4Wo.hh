// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Wo.hh,v 2.1 1998/11/06 13:42:01 allison Exp $
// GEANT4 tag $Name: geant4-00 $
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
  ~G4Wo ();
  G4VScene*    CreateScene (const G4String& name = "");
  G4VView*     CreateView  (G4VScene&, const G4String& name = "");
  static G4VInteractorManager* GetInteractorManager ();
};

#endif

#endif
