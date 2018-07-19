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
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// 
// G4UTrap
//
// Class description:
//
//   Wrapper class for G4Trap to make use of VecGeom Trapezoid.

// History:
// 13.09.13 G.Cosmo, CERN/PH
// --------------------------------------------------------------------
#ifndef G4UTrap_HH
#define G4UTrap_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <volumes/UnplacedTrapezoid.h>

#include "G4Polyhedron.hh"

class G4UTrap : public G4UAdapter<vecgeom::UnplacedTrapezoid> 
{
  using Shape_t = vecgeom::UnplacedTrapezoid;
  using Base_t = G4UAdapter<vecgeom::UnplacedTrapezoid>;

  public:  // with description

    G4UTrap( const G4String& pName,
                   G4double pDz,
                   G4double pTheta, G4double pPhi,
                   G4double pDy1, G4double pDx1, G4double pDx2,
                   G4double pAlp1,
                   G4double pDy2, G4double pDx3, G4double pDx4,
                   G4double pAlp2 );
      //
      // The most general constructor for G4Trap which prepares plane
      // equations and corner coordinates from parameters

    G4UTrap( const G4String& pName,
             const G4ThreeVector pt[8] ) ;
      //
      // Prepares plane equations and parameters from corner coordinates

    G4UTrap( const G4String& pName,
                   G4double pZ,
                   G4double pY,
                   G4double pX, G4double pLTX );
      //
      // Constructor for Right Angular Wedge from STEP (assumes pLTX<=pX)

    G4UTrap( const G4String& pName,
                   G4double pDx1,  G4double pDx2,
                   G4double pDy1,  G4double pDy2,
                   G4double pDz );
      //
      // Constructor for G4Trd       

    G4UTrap(const G4String& pName,
                  G4double pDx, G4double pDy, G4double pDz,
                  G4double pAlpha, G4double pTheta, G4double pPhi );
      //
      // Constructor for G4Para

    G4UTrap( const G4String& pName );
      //
      // Constructor for "nominal" G4Trap whose parameters are to be set
      // by a G4VPVParamaterisation later

   ~G4UTrap();

    void ComputeDimensions(      G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    G4VSolid* Clone() const;

    using Base_t::GetTanAlpha1;
    using Base_t::GetTanAlpha2;

    G4double GetZHalfLength()  const;
    G4double GetYHalfLength1() const;
    G4double GetXHalfLength1() const;
    G4double GetXHalfLength2() const;
    G4double GetYHalfLength2() const;
    G4double GetXHalfLength3() const;
    G4double GetXHalfLength4() const;
    G4double GetThetaCphi()    const;
    G4double GetThetaSphi()    const;
    TrapSidePlane GetSidePlane(G4int n) const;
    G4ThreeVector GetSymAxis() const;

    void SetAllParameters(G4double pDz, G4double pTheta, G4double pPhi,
                          G4double pDy1, G4double pDx1, G4double pDx2,
                          G4double pAlp1,
                          G4double pDy2, G4double pDx3, G4double pDx4,
                          G4double pAlp2);
    void SetPlanes(const G4ThreeVector pt[8]);

    inline G4GeometryType GetEntityType() const;

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                           G4double& pMin, G4double& pMax) const;

    G4Polyhedron* CreatePolyhedron() const;

  public:  // without description

    G4UTrap(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UTrap(const G4UTrap& rhs);
    G4UTrap& operator=(const G4UTrap& rhs); 
      // Copy constructor and assignment operator.
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UTrap::GetEntityType() const
{
  return "G4Trap";
}

#endif  // G4GEOM_USE_USOLIDS

#endif
