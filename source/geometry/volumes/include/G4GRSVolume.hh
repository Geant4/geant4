// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GRSVolume.hh,v 1.2 1999-12-15 14:50:21 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4GRSVolume     Paul Kent August 1996
//
// Object representing a touchable detector element - maintains
// associations between a physical volume and its net resultant
// local->global transform
//
// NOTE: The (optional) rotation matrix is copied

#ifndef G4GRSVOLUME_HH
#define G4GRSVOLUME_HH

#include "G4VTouchable.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"

class G4GRSVolume : public G4VTouchable
{
public:
  G4GRSVolume(G4VPhysicalVolume *pVol,
	      const G4RotationMatrix *pRot,
	      const G4ThreeVector &tlate);
  G4GRSVolume(G4VPhysicalVolume *pVol,
	      const G4RotationMatrix &rot,
	      const G4ThreeVector &tlate);
  ~G4GRSVolume();
  G4VPhysicalVolume* GetVolume(G4int depth=0) const;
  G4VSolid* GetSolid(G4int depth=0) const;
  const G4ThreeVector& GetTranslation(G4int depth=0) const;
  const G4RotationMatrix*  GetRotation(G4int depth=0) const;
  
private:
  
  G4VPhysicalVolume *fvol;
  G4RotationMatrix *frot;
  G4ThreeVector ftlate;
};

#include "G4GRSVolume.icc"
#endif
