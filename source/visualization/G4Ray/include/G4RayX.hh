// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RayX.hh,v 2.1 1998/11/06 13:41:52 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// Nikos Savvas  1st June 1997
// GEANT4 Ray Tracing graphics system factory.

#if defined (G4VIS_BUILD_RAYX_DRIVER) || defined (G4VIS_USE_RAYX)

#ifndef G4RAYX_HH
#define G4RAYX_HH

#include "G4VGraphicsSystem.hh"

class G4RayXMessenger;

class G4RayX: public G4VGraphicsSystem {
public:
  G4RayX ();
  ~G4RayX ();
  G4VScene* CreateScene (const G4String& name = "");
  G4VView*  CreateView  (G4VScene&, const G4String& name = "");
private:
  static G4int fInstances;
  G4RayXMessenger*   fpMessenger;        // Pointer to messenger.
};

#endif

#endif
