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
// G4UGenericPolycone
//
// Class description:
//
// Wrapper class for G4GenericPolycone to make use of VecGeom GenericPolycone.

// 30.10.13 G.Cosmo, CERN
// --------------------------------------------------------------------
#ifndef G4UGENERICPOLYCONE_HH
#define G4UGENERICPOLYCONE_HH

#include "G4UAdapter.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include <volumes/UnplacedGenericPolycone.h>

#include "G4TwoVector.hh"
#include "G4PolyconeSide.hh"

class G4UGenericPolycone : public G4UAdapter<vecgeom::UnplacedGenericPolycone>
{
  using Shape_t = vecgeom::UnplacedGenericPolycone;
  using Base_t  = G4UAdapter<vecgeom::UnplacedGenericPolycone>;

  public:  // with description

    G4UGenericPolycone(const G4String& name, 
                       G4double phiStart,    // initial phi starting angle
                       G4double phiTotal,    // total phi angle
                       G4int    numRZ,       // number corners in r,z space
                 const G4double r[],         // r coordinate of these corners
                 const G4double z[]       ); // z coordinate of these corners

   ~G4UGenericPolycone();

    G4double GetStartPhi()    const;
    G4double GetEndPhi()      const;
    G4double GetSinStartPhi() const;
    G4double GetCosStartPhi() const;
    G4double GetSinEndPhi()   const;
    G4double GetCosEndPhi()   const;
    G4bool IsOpen()           const;
    G4int  GetNumRZCorner()   const;
    G4PolyconeSideRZ GetCorner(G4int index) const;
  
    inline G4GeometryType GetEntityType() const;

    G4VSolid* Clone() const;

  public:  // without description

    G4UGenericPolycone(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UGenericPolycone( const G4UGenericPolycone& source );
    G4UGenericPolycone& operator=(const G4UGenericPolycone& source);
      // Copy constructor and assignment operator.

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                           G4double& pMin, G4double& pMax) const;

    G4Polyhedron* CreatePolyhedron() const;

  private:

    G4double wrStart;
    G4double wrDelta;
    std::vector<G4TwoVector> rzcorners;
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline G4GeometryType G4UGenericPolycone::GetEntityType() const
{
  return "G4GenericPolycone";
}
  
#endif  // G4GEOM_USE_USOLIDS

#endif
