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
//
// $Id: G4Tet.hh 104316 2017-05-24 13:04:23Z gcosmo $
//
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4Tet
//
// Class description:
//
//   A G4Tet is a tetrahedrasolid.
//

// History:
// -------
// 03.09.2004 - M.H.Mendenhall & R.A.Weller (Vanderbilt University, USA)
// 10.02.2005 - D.Anninos (CERN) - Added GetPointOnSurface() method.
// 12.11.2006 - M.H.Mendenhall - Added GetSurfaceArea() concrete implementation.
// 20.09.2010 - G.Cosmo (CERN) - Added copy-ctor and operator=().
// 24.09.2016 - E.Tcherniaev - Removed CreateRotatedVertices().
// --------------------------------------------------------------------
#ifndef G4TET_HH
#define G4TET_HH

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

  public:  // with description

    G4Tet(const G4String& pName, 
                G4ThreeVector anchor,
                G4ThreeVector p2,
                G4ThreeVector p3,
                G4ThreeVector p4, 
                G4bool *degeneracyFlag=0);

    virtual ~G4Tet();

    void ComputeDimensions(G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pmin, G4double& pmax) const;
    // Methods for solid

    EInside Inside(const G4ThreeVector& p) const;

    G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const;

    G4double DistanceToIn(const G4ThreeVector& p, const G4ThreeVector& v) const;

    G4double DistanceToIn(const G4ThreeVector& p) const;

    G4double DistanceToOut(const G4ThreeVector& p, const G4ThreeVector& v,
                           const G4bool calcNorm=false,
                                 G4bool *validNorm=0, G4ThreeVector *n=0) const;

    G4double DistanceToOut(const G4ThreeVector& p) const;

    G4double GetCubicVolume();
    G4double GetSurfaceArea();

    G4GeometryType GetEntityType() const;

    G4VSolid* Clone() const;

    std::ostream& StreamInfo(std::ostream& os) const;

    G4ThreeVector GetPointOnSurface() const;

    // Functions for visualization

    void          DescribeYourselfTo (G4VGraphicsScene& scene) const;
    G4VisExtent   GetExtent          () const;
    G4Polyhedron* CreatePolyhedron   () const;
    G4Polyhedron* GetPolyhedron      () const;

  public:   // without description

    G4Tet(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4Tet(const G4Tet& rhs);
    G4Tet& operator=(const G4Tet& rhs); 
      // Copy constructor and assignment operator.

    const char* CVSHeaderVers()
      { return "$Id: G4Tet.hh 104316 2017-05-24 13:04:23Z gcosmo $"; }
    const char* CVSFileVers()
      { return CVSVers; }
    void PrintWarnings(G4bool flag)
      { warningFlag=flag; }
    static G4bool CheckDegeneracy(G4ThreeVector anchor,
                                  G4ThreeVector p2,
                                  G4ThreeVector p3,
                                  G4ThreeVector p4);
    std::vector<G4ThreeVector> GetVertices() const;
      // Return the four vertices of the shape.

  private:

    G4double fCubicVolume, fSurfaceArea;

    mutable G4bool fRebuildPolyhedron;
    mutable G4Polyhedron* fpPolyhedron;

    G4ThreeVector GetPointOnFace(G4ThreeVector p1, G4ThreeVector p2, 
                                 G4ThreeVector p3, G4double& area) const;
    static const char CVSVers[];

  private:

    G4ThreeVector fAnchor, fP2, fP3, fP4, fMiddle;
    G4ThreeVector fNormal123, fNormal142, fNormal134, fNormal234;

    G4bool warningFlag;

    G4double fCdotN123, fCdotN142, fCdotN134, fCdotN234;
    G4double fXMin, fXMax, fYMin, fYMax, fZMin, fZMax;
    G4double fDx, fDy, fDz, fTol, fMaxSize;
};

#endif

#endif
