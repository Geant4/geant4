//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4GRSVolume.hh,v 1.5 2001-07-11 10:00:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4GRSVolume
//
// Class description:
//
// Object representing a touchable detector element - maintains
// associations between a physical volume and its net resultant
// local->global transform.
//
// NOTE: The (optional) rotation matrix is copied

// History:
// - Created. Paul Kent, August 1996

#ifndef G4GRSVOLUME_HH
#define G4GRSVOLUME_HH

#include "G4VTouchable.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"

class G4GRSVolume : public G4VTouchable
{
  public:  // with description

    G4GRSVolume(G4VPhysicalVolume *pVol,
	        const G4RotationMatrix *pRot,
	        const G4ThreeVector &tlate);
    G4GRSVolume(G4VPhysicalVolume *pVol,
	        const G4RotationMatrix &rot,
	        const G4ThreeVector &tlate);
    ~G4GRSVolume();

    inline G4VPhysicalVolume* GetVolume(G4int depth=0) const;
    inline G4VSolid* GetSolid(G4int depth=0) const;
    inline const G4ThreeVector& GetTranslation(G4int depth=0) const;
    inline const G4RotationMatrix*  GetRotation(G4int depth=0) const;

  private:

    G4GRSVolume(const G4GRSVolume&);
    G4GRSVolume& operator=(const G4GRSVolume&);
      // Copy constructor and assignment operator NOT public

  private:
  
    G4VPhysicalVolume *fvol;
    G4RotationMatrix *frot;
    G4ThreeVector ftlate;
};

#include "G4GRSVolume.icc"

#endif
