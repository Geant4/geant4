// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GRSSolid.hh,v 1.3 2000-04-25 16:15:02 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4GRSSolid
//
// Class description:
//
// Object representing a touchable solid - maintains the association
// between a solid and its net resultant local->global transform.
//
// NOTE: The (optional) rotation matrix is copied

// History:
// - Created. Paul Kent, August 1996

#ifndef G4GRSSOLID_HH
#define G4GRSSOLID_HH

#include "G4VTouchable.hh"

#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"

class G4GRSSolid : public G4VTouchable
{
  public:  // with description

    G4GRSSolid(G4VSolid *pSolid,
	       const G4RotationMatrix *pRot,
	       const G4ThreeVector &tlate);
    G4GRSSolid(G4VSolid *pSolid,
	       const G4RotationMatrix &rot,
	       const G4ThreeVector &tlate);
    ~G4GRSSolid();

    G4VSolid* GetSolid(G4int depth=0) const;
    const G4ThreeVector& GetTranslation(G4int depth=0) const;
    const G4RotationMatrix*  GetRotation(G4int depth=0) const;
  
  private:
  
    G4VSolid *fsolid;
    G4RotationMatrix *frot;
    G4ThreeVector ftlate;
};

#include "G4GRSSolid.icc"

#endif
