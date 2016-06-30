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
// class G4BoundingEnvelope
//
// Class description:
//
// Helper class to facilitate calculation of the extent of a solid
// within the limits defined by the G4VoxelLimits object.
//
// The function CalculateExtent() of a particular solid can create
// a G4BoundingEnvelope object that bounds the solid and then call
// CalculateExtent() of the G4BoundingEnvelope object.
//
// Calculation of the extent by G4BoundingEnvelope takes into account
// special parameter "delta" - width of an imagery layer that envelops
// the object. This value, multiplied by max scale factor, is added
// to the voxel limits during calculation of the extent.
//
// Example of use. 
// In case of G4Box, max possible distance of a point to the border
// of the box, where the point is still considered as belonging to 
// the surface of the box, is sqrt(0.5)*kCarTolerance (see corner).
// So, it will be safe to set the extension = kCarTolerance.
//
// Alternative solution can be to define G4BoundingEnvelope wide
// enough to include the surface of the solid and set delta = 0. 
//
// The class supports the following bouding envelopes:
// - axis aligned bounding box (AABB);
// - bounding prism given by two convex polygonal bases;
// - bounding pyramid given by apex and convex polygonal base;
// - set of bounding prisms, given by a sequence of convex
//   polygonal bases;

// History:
//
// 2016.05.25 E.Tcherniaev - initial version
//
// --------------------------------------------------------------------
#ifndef G4BOUNDINGENVELOPE_HH
#define G4BOUNDINGENVELOPE_HH

#include <vector>
#include "geomdefs.hh"

#include "G4ThreeVector.hh"
#include "G4VoxelLimits.hh"
#include "G4Transform3D.hh"
#include "G4Point3D.hh"
#include "G4Plane3D.hh"

typedef std::vector<G4ThreeVector>     G4ThreeVectorList;
typedef std::vector<G4Point3D>         G4Polygon3D;
typedef std::pair<G4Point3D,G4Point3D> G4Segment3D;

class G4BoundingEnvelope
{
  public:

    G4BoundingEnvelope(const G4ThreeVector& pMin, 
                       const G4ThreeVector& pMax, G4double delta);
      // Constructor from an axis aligned bounding box (AABB)

    G4BoundingEnvelope(const G4ThreeVectorList& baseA,
                       const G4ThreeVectorList& baseB, G4double delta);
      // Constructor from a prism given by two bases, the bases
      // should have equal number of vertices

    G4BoundingEnvelope(const G4ThreeVector& apex,
                       const G4ThreeVectorList& base, G4double delta);
      // Constructor from a pyramid given by apex and base

    G4BoundingEnvelope(const std::vector<G4ThreeVectorList*>& polygons,
                       G4double delta);
      // Constructor from a sequence of convex polygons, the polygons
      // should have equal numbers of vertices except first and last
      // polygons which may consist of a single vertex

    G4BoundingEnvelope(const G4BoundingEnvelope& rhs);
      // Copy constructor

    G4BoundingEnvelope& operator=(const G4BoundingEnvelope& rhs);
      // Assignment operator

    ~G4BoundingEnvelope();
      // Destructor

    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimits,
                           const G4Transform3D& pTransform3D,
                           G4double& pMin, G4double& pMax) const;
      // Calculate extent of the envelope

  private:

    void SetDelta(G4double delta);
      // Set the extension

    void SetBoundingBox(const G4ThreeVector& pMin,
                        const G4ThreeVector& pMax);
      // Set AABB (axis aligned bounding box)

    void SetBoundingPrism(const G4ThreeVectorList& baseA,
                          const G4ThreeVectorList& baseB);
      // Set bounding prism

    void SetBoundingPyramid(const G4ThreeVector& apex,
                            const G4ThreeVectorList& base);
      // Set bounding pyramid

    void SetBoundingPolygons(const std::vector<G4ThreeVectorList*>& polygons);
      // Set bounding sequence of convex polygons

    void CleanPolygons();
      // Free allocated memory

    G4VoxelLimits GetAdjustedVoxelLimits(const G4VoxelLimits& pVoxelLimits,
                                         G4double pDelta) const;
      // Extend voxel limits by scaled surface tolerance 

    void TransformVertices(const G4Transform3D& pTransform3D,
                           const G4Polygon3D& polyA,
                                 G4Polygon3D& polyB,
                                 G4Segment3D& pAABB) const;
      // Transform vertices of a polygon and update AABB (bounding box)

    void CreateListOfEdges(const G4Polygon3D& baseA,
                           const G4Polygon3D& baseB,
                           std::vector<G4Segment3D>& pEdges) const;
      // Create list of edges of a prism
  
    void CreateListOfPlanes(const G4Polygon3D& baseA,
                            const G4Polygon3D& baseB,
                            std::vector<G4Plane3D>& pPlanes) const;
      // Create list of planes bounding a prism

    G4bool ClipEdgesByVoxelLimits(const std::vector<G4Segment3D>& pEdges,
                                  const G4VoxelLimits& pLimits,
                                        G4Segment3D& pExtent) const;
      // Clip set of edges by G4VoxelLimits

    void ClipVoxelLimitsByPlanes(const G4VoxelLimits& pLimits,
                                 const std::vector<G4Plane3D>& pPlanes,
                                 const G4Segment3D& pAABB,
                                       G4Segment3D& pExtent) const;
      // Clip G4VoxelLimits by set of planes bounding a convex prism

  private:

    G4double fDelta;                  // extention
    std::vector<G4Polygon3D*> fBases; // sequence of polygonal bases
};

#endif // G4BOUNDINGENVELOPE_HH
