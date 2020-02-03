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
// class G4GRSSolid
//
// Class description:
//
// Object representing a touchable solid - maintains the association
// between a solid and its net resultant local->global transform.
//
// NOTE: The (optional) rotation matrix is copied

// Created: Paul Kent - August 1996
// ----------------------------------------------------------------------
#ifndef G4GRSSOLID_HH
#define G4GRSSOLID_HH

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VTouchable.hh"

class G4VSolid;

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

    G4GRSSolid(const G4GRSSolid&) = delete;
    G4GRSSolid& operator=(const G4GRSSolid&) = delete;
      // Copy constructor and assignment operator not allowed

    inline G4VSolid* GetSolid(G4int depth=0) const;
    inline const G4ThreeVector& GetTranslation(G4int depth=0) const;
    inline const G4RotationMatrix*  GetRotation(G4int depth=0) const;

  private:
  
    G4VSolid* fsolid = nullptr;
    G4RotationMatrix* frot = nullptr;
    G4ThreeVector ftlate;
};

#include "G4GRSSolid.icc"

#endif
