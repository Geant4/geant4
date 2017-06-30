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
// * technical work of the GEANT4 collaboration and of QinetiQ Ltd,   *
// * subject to DEFCON 705 IPR conditions.                            *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4TessellatedSolid.hh 104316 2017-05-24 13:04:23Z gcosmo $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// Class G4TessellatedSolid
//
// Class description:
//
//    G4TessellatedSolid is a special Geant4 solid defined by a number of 
//    facets (UVFacet). It is important that the supplied facets shall form a
//    fully enclose space which is the solid. 
//    At the moment only two types of facet can be used for the construction of 
//    a G4TessellatedSolid, i.e. the G4TriangularFacet and G4QuadrangularFacet.
//
//    How to contruct a G4TessellatedSolid:
//  
//    First declare a tessellated solid:
//
//      G4TessellatedSolid* solidTarget = new G4TessellatedSolid("Solid_name");
//
//    Define the facets which form the solid
// 
//      G4double targetSiz = 10*cm ;
//      G4TriangularFacet *facet1 = new
//      G4TriangularFacet (G4ThreeVector(-targetSize,-targetSize,        0.0),
//                         G4ThreeVector(+targetSize,-targetSize,        0.0),
//                         G4ThreeVector(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      G4TriangularFacet *facet2 = new
//      G4TriangularFacet (G4ThreeVector(+targetSize,-targetSize,        0.0),
//                         G4ThreeVector(+targetSize,+targetSize,        0.0),
//                         G4ThreeVector(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      G4TriangularFacet *facet3 = new
//      G4TriangularFacet (G4ThreeVector(+targetSize,+targetSize,        0.0),
//                         G4ThreeVector(-targetSize,+targetSize,        0.0),
//                         G4ThreeVector(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      G4TriangularFacet *facet4 = new
//      G4TriangularFacet (G4ThreeVector(-targetSize,+targetSize,        0.0),
//                         G4ThreeVector(-targetSize,-targetSize,        0.0),
//                         G4ThreeVector(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      G4QuadrangularFacet *facet5 = new
//      G4QuadrangularFacet (G4ThreeVector(-targetSize,-targetSize,      0.0),
//                           G4ThreeVector(-targetSize,+targetSize,      0.0),
//                           G4ThreeVector(+targetSize,+targetSize,      0.0),
//                           G4ThreeVector(+targetSize,-targetSize,      0.0),
//                           ABSOLUTE);
//
//    Then add the facets to the solid:    
//
//      solidTarget->AddFacet((UVFacet*) facet1);
//      solidTarget->AddFacet((UVFacet*) facet2);
//      solidTarget->AddFacet((UVFacet*) facet3);
//      solidTarget->AddFacet((UVFacet*) facet4);
//      solidTarget->AddFacet((UVFacet*) facet5);
//
//    Finally declare the solid is complete:
//
//      solidTarget->SetSolidClosed(true);

// CHANGE HISTORY
// --------------
// 31 October 2004, P R Truscott, QinetiQ Ltd, UK
//  - Created.
// 22 November 2005, F Lei, 
//  - Added GetPolyhedron().
// 12 October 2012, M Gayer,
//  - Reviewed optimized implementation including voxelization of surfaces.
//
///////////////////////////////////////////////////////////////////////////////
#ifndef G4TessellatedSolid_hh
#define G4TessellatedSolid_hh 1

#include <iostream>
#include <vector>
#include <set>
#include <map>

#include "G4VSolid.hh"
#include "G4Types.hh"
#include "G4Voxelizer.hh"

struct G4VertexInfo
{
  G4int id;
  G4double mag2;
};

class G4VFacet;

class G4VertexComparator
{
public:
  G4bool operator() (const G4VertexInfo &l, const G4VertexInfo &r) const
  {
    return l.mag2 == r.mag2 ? l.id < r.id : l.mag2 < r.mag2;
  }
};

class G4TessellatedSolid : public G4VSolid
{
  public:  // with description

    G4TessellatedSolid ();
    virtual ~G4TessellatedSolid ();

    G4TessellatedSolid (const G4String &name);

    G4TessellatedSolid(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4TessellatedSolid (const G4TessellatedSolid &ts);
    G4TessellatedSolid &operator= (const G4TessellatedSolid &right);
    G4TessellatedSolid &operator+= (const G4TessellatedSolid &right);

    G4bool AddFacet (G4VFacet *aFacet);
    inline G4VFacet *GetFacet (G4int i) const;

    G4int GetNumberOfFacets () const;

    virtual EInside Inside (const G4ThreeVector &p) const;
    virtual G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const;
    virtual G4double DistanceToIn(const G4ThreeVector& p,
                                  const G4ThreeVector& v)const;
    virtual G4double DistanceToIn(const G4ThreeVector& p) const;
    virtual G4double DistanceToOut(const G4ThreeVector& p) const;
    virtual G4double DistanceToOut(const G4ThreeVector& p,
                                   const G4ThreeVector& v,
                                   const G4bool calcNorm,
                                         G4bool *validNorm,
                                         G4ThreeVector *norm) const;

    virtual G4bool Normal (const G4ThreeVector &p, G4ThreeVector &n) const;
    virtual G4double SafetyFromOutside(const G4ThreeVector &p,
                                             G4bool aAccurate=false) const;
    virtual G4double SafetyFromInside (const G4ThreeVector &p,
                                             G4bool aAccurate=false) const;

    virtual G4GeometryType GetEntityType () const;
    virtual std::ostream &StreamInfo(std::ostream &os) const;

    virtual G4VSolid* Clone() const;

    virtual G4ThreeVector GetPointOnSurface() const;
    virtual G4double GetSurfaceArea();
    virtual G4double GetCubicVolume ();

    void SetSolidClosed (const G4bool t);
    G4bool GetSolidClosed () const;

    inline void SetMaxVoxels(G4int max);

    inline G4Voxelizer &GetVoxels();

    virtual G4bool CalculateExtent(const EAxis pAxis,
                                   const G4VoxelLimits& pVoxelLimit,
                                   const G4AffineTransform& pTransform,
                                         G4double& pMin, G4double& pMax) const;

    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;

    G4double      GetMinXExtent () const;
    G4double      GetMaxXExtent () const;
    G4double      GetMinYExtent () const;
    G4double      GetMaxYExtent () const;
    G4double      GetMinZExtent () const;
    G4double      GetMaxZExtent () const;

    virtual G4Polyhedron* CreatePolyhedron () const;
    virtual G4Polyhedron* GetPolyhedron    () const;
    virtual void DescribeYourselfTo (G4VGraphicsScene& scene) const;
    virtual G4VisExtent   GetExtent () const;

    G4int AllocatedMemoryWithoutVoxels();
    G4int AllocatedMemory();
    void DisplayAllocatedMemory();

  private: // without description

    void Initialize();

    G4double DistanceToOutNoVoxels(const G4ThreeVector &p,
                                   const G4ThreeVector &v,
                                         G4ThreeVector &aNormalVector,
                                         G4bool        &aConvex,
                                         G4double aPstep = kInfinity) const;
    G4double DistanceToInCandidates(const std::vector<G4int> &candidates,
                                    const G4ThreeVector &aPoint,
                                    const G4ThreeVector &aDirection) const;
    void DistanceToOutCandidates(const std::vector<G4int> &candidates,
                                 const G4ThreeVector &aPoint,
                                 const G4ThreeVector &direction,
                                       G4double &minDist,
                                       G4ThreeVector &minNormal,
                                       G4int &minCandidate) const;
    G4double DistanceToInNoVoxels(const G4ThreeVector &p,
                                  const G4ThreeVector &v,
                                        G4double aPstep = kInfinity) const;
    void SetExtremeFacets();

    EInside InsideNoVoxels (const G4ThreeVector &p) const;
    EInside InsideVoxels(const G4ThreeVector &aPoint) const;

    void Voxelize();

    void CreateVertexList();

    void PrecalculateInsides();

    void SetRandomVectors();

    G4double DistanceToInCore(const G4ThreeVector &p, const G4ThreeVector &v,
                                    G4double aPstep = kInfinity) const;
    G4double DistanceToOutCore(const G4ThreeVector &p, const G4ThreeVector &v,
                                     G4ThreeVector &aNormalVector,
                                     G4bool        &aConvex,
                                     G4double aPstep = kInfinity) const;

    G4int SetAllUsingStack(const std::vector<G4int> &voxel,
                           const std::vector<G4int> &max,
                                 G4bool status, G4SurfBits &checked);

    void DeleteObjects ();
    void CopyObjects (const G4TessellatedSolid &s);

    static G4bool CompareSortedVoxel(const std::pair<G4int, G4double> &l,
                                     const std::pair<G4int, G4double> &r);

    G4double MinDistanceFacet(const G4ThreeVector &p, G4bool simple,
                                    G4VFacet * &facet) const;

    inline G4bool OutsideOfExtent(const G4ThreeVector &p,
                                        G4double tolerance=0) const;

  protected:

    G4double kCarToleranceHalf;

  private:

    mutable G4bool fRebuildPolyhedron;
    mutable G4Polyhedron* fpPolyhedron;

    std::vector<G4VFacet *>  fFacets;
    std::set<G4VFacet *> fExtremeFacets; // Does all other facets lie on
                                         // or behind this surface?

    G4GeometryType           fGeometryType;
    G4double                 fCubicVolume;
    G4double                 fSurfaceArea;

    std::vector<G4ThreeVector>  fVertexList;

    std::set<G4VertexInfo,G4VertexComparator> fFacetList;

    G4ThreeVector fMinExtent, fMaxExtent;

    G4bool fSolidClosed;

    std::vector<G4ThreeVector> fRandir;

    G4int fMaxTries;

    G4Voxelizer fVoxels;  // Pointer to the voxelized solid

    G4SurfBits fInsides;
};

///////////////////////////////////////////////////////////////////////////////
// Inlined Methods
///////////////////////////////////////////////////////////////////////////////

inline G4VFacet *G4TessellatedSolid::GetFacet (G4int i) const
{
  return fFacets[i];
}

inline void G4TessellatedSolid::SetMaxVoxels(G4int max)
{
  fVoxels.SetMaxVoxels(max);
}

inline G4Voxelizer &G4TessellatedSolid::GetVoxels()
{
  return fVoxels;
}

inline G4bool G4TessellatedSolid::OutsideOfExtent(const G4ThreeVector &p,
                                                  G4double tolerance) const
{
  return ( p.x() < fMinExtent.x() - tolerance
        || p.x() > fMaxExtent.x() + tolerance
        || p.y() < fMinExtent.y() - tolerance
        || p.y() > fMaxExtent.y() + tolerance
        || p.z() < fMinExtent.z() - tolerance
        || p.z() > fMaxExtent.z() + tolerance);
}

#endif
