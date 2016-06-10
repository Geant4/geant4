//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4GRSVolume.hh 66356 2012-12-18 09:02:32Z gcosmo $
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
// ----------------------------------------------------------------------
#ifndef G4GRSVOLUME_HH
#define G4GRSVOLUME_HH

#include "G4VTouchable.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
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
