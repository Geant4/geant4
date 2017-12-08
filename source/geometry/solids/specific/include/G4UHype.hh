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
// $Id:$
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4UHype
//
// Class description:
//
//   Wrapper class for G4Hype to make use of VecGeom Hyperboloid.

// History:
// 16.10.17 G.Cosmo, CERN
// --------------------------------------------------------------------
#ifndef G4UHYPE_HH
#define G4UHYPE_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <volumes/UnplacedHype.h>

#include "G4Polyhedron.hh"

class G4UHype : public G4UAdapter<vecgeom::UnplacedHype>
{
  using Shape_t = vecgeom::UnplacedHype;
  using Base_t  = G4UAdapter<vecgeom::UnplacedHype>;

  public:  // with description

    G4UHype(const G4String& name,
                  G4double  newInnerRadius,
                  G4double  newOuterRadius,
                  G4double  newInnerStereo,
                  G4double  newOuterStereo,
                  G4double  newHalfLenZ);
   ~G4UHype();

    void ComputeDimensions(      G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    G4VSolid* Clone() const;

    G4double GetInnerRadius () const;
    G4double GetOuterRadius () const;
    G4double GetZHalfLength () const;
    G4double GetInnerStereo () const;
    G4double GetOuterStereo () const;

    void SetInnerRadius (G4double newIRad);
    void SetOuterRadius (G4double newORad);
    void SetZHalfLength (G4double newHLZ);
    void SetInnerStereo (G4double newISte);
    void SetOuterStereo (G4double newOSte);

    inline G4GeometryType GetEntityType() const;

  public:  // without description

    G4UHype(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UHype( const G4UHype &source );
    G4UHype &operator=( const G4UHype &source );
      // Copy constructor and assignment operator.

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;
    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                           G4double& pmin, G4double& pmax) const;
    G4Polyhedron* CreatePolyhedron() const;
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UHype::GetEntityType() const
{
  return "G4Hype";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
