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
// G4UPara
//
// Class description:
//
//   Wrapper class for G4Para to make use of VecGeom Parallelepiped.

// History:
// 13.09.13 G.Cosmo, CERN/PH
// --------------------------------------------------------------------

#ifndef G4UPARA_HH
#define G4UPARA_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <volumes/UnplacedParallelepiped.h>

#include "G4Polyhedron.hh"

class G4UPara : public G4UAdapter<vecgeom::UnplacedParallelepiped> 
{
  using Shape_t = vecgeom::UnplacedParallelepiped;
  using Base_t  = G4UAdapter<vecgeom::UnplacedParallelepiped>;

  public:  // with description

    G4UPara(const G4String& pName,
                  G4double pDx, G4double pDy, G4double pDz,
                  G4double pAlpha, G4double pTheta, G4double pPhi);

    G4UPara(const G4String& pName,
            const G4ThreeVector pt[8]);

   ~G4UPara();

    // Accessors

    G4double GetZHalfLength()  const;
    G4double GetYHalfLength()  const;
    G4double GetXHalfLength()  const;
    G4ThreeVector GetSymAxis() const;
    G4double GetTanAlpha()     const;

    // Modifiers

    void SetXHalfLength(G4double val);
    void SetYHalfLength(G4double val);
    void SetZHalfLength(G4double val);
    void SetAlpha(G4double alpha);
    void SetTanAlpha(G4double val);
    void SetThetaAndPhi(double pTheta, double pPhi);

    void SetAllParameters(G4double pDx, G4double pDy, G4double pDz,
                          G4double pAlpha, G4double pTheta, G4double pPhi);

    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    inline G4GeometryType GetEntityType() const;

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const;

    G4VSolid* Clone() const;

    G4Polyhedron* CreatePolyhedron   () const;

  public:  // without description

    G4UPara(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects

    G4UPara(const G4UPara& rhs);
    G4UPara& operator=(const G4UPara& rhs);
      // Copy constructor and assignment operator

  private:

    void CheckParameters();
      // Check parameters

    void MakePlanes();
      // Set side planes

  private:

    G4double fTalpha,fTthetaCphi,fTthetaSphi;
    struct { G4double a,b,c,d; } fPlanes[4];
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UPara::GetEntityType() const
{
  return "G4Para";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
