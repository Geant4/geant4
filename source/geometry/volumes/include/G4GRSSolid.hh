// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GRSSolid.hh,v 1.1 1999-01-07 16:08:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4GRSSolid     Paul Kent August 1996
//
// Object representing a touchable solid - maintains
// association between a solid and its net resultant
// local->global transform
//
// NOTE: The (optional) rotation matrix is copied

#ifndef G4GRSSOLID_HH
#define G4GRSSOLID_HH

#include "G4VTouchable.hh"

#include "G4VPhysicalVolume.hh"
// #include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"

class G4GRSSolid : public G4VTouchable
{
public:
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
