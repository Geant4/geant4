// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Xo.hh,v 2.1 1998/11/06 13:42:04 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// Guy Barrand 04 November 1996
// Wo graphics system factory.

#ifndef G4XO_HH
#define G4XO_HH

#if defined(G4VIS_BUILD_OPACS_DRIVER) || defined(G4VIS_USE_OPACS)

//G4
#include "G4VGraphicsSystem.hh"

class G4VInteractorManager;

class G4Xo: public G4VGraphicsSystem {
public:
  G4Xo ();                     
  ~G4Xo ();                     
  G4VScene*           CreateScene   (const G4String& name = "");
  G4VView*            CreateView    (G4VScene&, const G4String& name = "");
  static G4VInteractorManager* GetInteractorManager ();
};

#endif

#endif
