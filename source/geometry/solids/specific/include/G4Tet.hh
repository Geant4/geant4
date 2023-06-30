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
// * This  code  implementation is the  intellectual property  of the *
// * Vanderbilt University Free Electron Laser Center                 *
// * Vanderbilt University, Nashville, TN, USA                        *
// * Development supported by:                                        *
// * United States MFEL program  under grant FA9550-04-1-0045         *
// * and NASA under contract number NNG04CT05P.                       *
// * Written by Marcus H. Mendenhall and Robert A. Weller.            *
// *                                                                  *
// * Contributed to the Geant4 Core, January, 2005.                   *
// *                                                                  *
// ********************************************************************
//
// G4Tet
//
// Class description:
//
// A G4Tet is a tetrahedra solid.

// 03.09.2004 - M.H.Mendenhall & R.A.Weller (Vanderbilt University, USA)
// 08.01.2020 - E.Tcherniaev, complete revision, speed up
// --------------------------------------------------------------------
#ifndef G4TET_HH
#define G4TET_HH

#include "G4GeomTypes.hh"

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UTET 1
#endif

#if defined(G4GEOM_USE_UTET)
  #define G4UTet G4Tet
  #include "G4UTet.hh"
#else

#include "G4VSolid.hh"

class G4Tet : public G4VSolid
{
  public:

    G4Tet(const G4String& pName,
          const G4ThreeVector& anchor,
          const G4ThreeVector& p1,
          const G4ThreeVector& p2,
          const G4ThreeVector& p3,
                G4bool* degeneracyFlag = nullptr);

    ~G4Tet() override;

    void SetVertices(const G4ThreeVector& anchor,
                     const G4ThreeVector& p1,
                     const G4ThreeVector& p2,
                     const G4ThreeVector& p3,
                     G4bool* degeneracyFlag = nullptr);

    // Accessors, return the four vertices of the shape
    void GetVertices(G4ThreeVector& anchor,
                     G4ThreeVector& p1,
                     G4ThreeVector& p2,
                     G4ThreeVector& p3) const;
    std::vector<G4ThreeVector> GetVertices() const;

    // Set warning flag - deprecated (dummy)
    inline void PrintWarnings(G4bool) {};

    // Return true if the tetrahedron is degenerate
    G4bool CheckDegeneracy(const G4ThreeVector& p0,
                           const G4ThreeVector& p1,
                           const G4ThreeVector& p2,
                           const G4ThreeVector& p3) const;

    // Standard methods
    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep) override;

    void SetBoundingLimits(const G4ThreeVector& pMin, const G4ThreeVector& pMax);
    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const override;
    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pmin, G4double& pmax) const override;

    EInside Inside(const G4ThreeVector& p) const override;
    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const override;
    G4double DistanceToIn(const G4ThreeVector& p,
                          const G4ThreeVector& v) const override;
    G4double DistanceToIn(const G4ThreeVector& p) const override;
    G4double DistanceToOut(const G4ThreeVector& p,
                           const G4ThreeVector& v,
                           const G4bool calcNorm = false,
                                 G4bool* validNorm = nullptr,
                                 G4ThreeVector* n = nullptr) const override;
    G4double DistanceToOut(const G4ThreeVector& p) const override;

    G4GeometryType GetEntityType() const override;

    G4VSolid* Clone() const override;

    std::ostream& StreamInfo(std::ostream& os) const override;

    G4double GetCubicVolume() override;
    G4double GetSurfaceArea() override;

    G4ThreeVector GetPointOnSurface() const override;

    // Methods for visualization
    void DescribeYourselfTo (G4VGraphicsScene& scene) const override;
    G4VisExtent GetExtent () const override;
    G4Polyhedron* CreatePolyhedron () const override;
    G4Polyhedron* GetPolyhedron () const override;

    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects
    G4Tet(__void__&);

    // Copy constructor
    G4Tet(const G4Tet& rhs);

    // Assignment operator
    G4Tet& operator=(const G4Tet& rhs);

  private:

    // Set data members
    void Initialize(const G4ThreeVector& p0,
                    const G4ThreeVector& p1,
                    const G4ThreeVector& p2,
                    const G4ThreeVector& p3);

    // Return normal to surface closest to p
    G4ThreeVector ApproxSurfaceNormal(const G4ThreeVector& p) const;

  private:

    G4double halfTolerance = 0;
    G4double fCubicVolume = 0; // Volume
    G4double fSurfaceArea = 0; // Surface area
    mutable G4bool fRebuildPolyhedron = false;
    mutable G4Polyhedron* fpPolyhedron = nullptr;

    G4ThreeVector fVertex[4];   // thetrahedron vertices
    G4ThreeVector fNormal[4];   // normals to faces
    G4double fDist[4] = {0};    // distances from origin to faces
    G4double fArea[4] = {0};    // face areas
    G4ThreeVector fBmin, fBmax; // bounding box
};

#endif

#endif
