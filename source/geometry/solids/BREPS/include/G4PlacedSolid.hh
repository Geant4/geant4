// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PlacedSolid.hh,v 1.1 1999-01-07 16:07:33 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4PLACEDSOLID_HH
#define G4PLACEDSOLID_HH

#include "G4BREPSolid.hh"


class G4PlacedSolid
{
public:

  G4PlacedSolid();
  G4PlacedSolid(G4BREPSolid*, G4Axis2Placement3D* =0);
  ~G4PlacedSolid();
  
  G4VSolid*      GetSolid()       { return solid;            }
  HepRotation*   GetRotation()    { return solidRotation;    }
  G4ThreeVector* GetTranslation() { return solidTranslation; } 

  G4bool operator==(const G4PlacedSolid& ps) const 
  {
    return (this==&ps) ? true : false; 
  }

private:

  G4BREPSolid*   solid;
  HepRotation*   solidRotation;
  G4ThreeVector* solidTranslation;

};

#endif




