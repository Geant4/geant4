// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VGraphicsSystem.hh,v 1.6 1999-12-15 14:54:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  27th March 1996
//
// Class description
//
// Abstract interface class for graphics systems.

#ifndef G4VGRAPHICSSYSTEM_HH
#define G4VGRAPHICSSYSTEM_HH

#include "globals.hh"

class G4VSceneHandler;
class G4VViewer;

class G4VGraphicsSystem {

public: // With description

  enum Functionality {
    noFunctionality,
    twoD,               // Simple 2D, e.g., X (no stored structures).
    twoDStore,          // 2D with stored structures.
    threeD,             // Passive 3D (with stored structures).
    threeDInteractive,  // 3D with "pick" functionality.
    virtualReality      // Virtual Reality functionality.
  };

  G4VGraphicsSystem (const G4String& name,
		     Functionality f);

  G4VGraphicsSystem (const G4String& name,
		     const G4String& nickname,
		     Functionality f);

  G4VGraphicsSystem (const G4String& name,
		     const G4String& nickname,
		     const G4String& description,
		     Functionality f);

  virtual ~G4VGraphicsSystem ();

  // For G4RWTPtrOrderedVector...
  G4bool operator == (const G4VGraphicsSystem& system) const;

  virtual G4VSceneHandler* CreateSceneHandler (const G4String& name) = 0;

  virtual G4VViewer* CreateViewer (G4VSceneHandler&, const G4String& name) = 0;

  // Access functions.
  const G4String& GetName          () const;
  const G4String& GetNickname      () const;
  const G4String& GetDescription   () const;
  Functionality   GetFunctionality () const;

protected:
  const G4String fName;
  const G4String fNickname;
  const G4String fDescription;
  Functionality  fFunctionality;
};

G4std::ostream& operator << (G4std::ostream& os, const G4VGraphicsSystem& gs);

#include "G4VGraphicsSystem.icc"

#endif
