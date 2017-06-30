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
// $Id: G4GenericTrap.hh 104316 2017-05-24 13:04:23Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4GenericTrap
//
// Class description:
//
// G4GenericTrap is a solid which represents an arbitrary trapezoid with 
// up to 8 vertices standing on two parallel planes perpendicular to Z axis.
// 
// Parameters in the constructor:
// - name               - solid name
// - halfZ              - the solid half length in Z
// - vertices           - the (x,y) coordinates of vertices:
//                        o first four points: vertices[i], i<4 
//                          are the vertices sitting on the -halfZ plane;
//                        o last four points: vertices[i], i>=4
//                          are the vertices sitting on the +halfZ plane.
//
//   The order of defining the vertices of the solid is the following:
//      - point 0 is connected with points 1,3,4
//      - point 1 is connected with points 0,2,5
//      - point 2 is connected with points 1,3,6
//      - point 3 is connected with points 0,2,7
//      - point 4 is connected with points 0,5,7
//      - point 5 is connected with points 1,4,6
//      - point 6 is connected with points 2,5,7
//      - point 7 is connected with points 3,4,6
// Points can be identical in order to create shapes with less than
// 8 vertices.

// Authors:
//   Tatiana Nikitina, CERN; Ivana Hrivnacova, IPN Orsay
//   Adapted from Root Arb8 implementation, author Andrea Gheata, CERN
// -------------------------------------------------------------------

#ifndef G4GenericTrap_HH
#define G4GenericTrap_HH

#if defined(G4GEOM_USE_USOLIDS)
#define G4GEOM_USE_UGENERICTRAP 1
#endif

#if defined(G4GEOM_USE_UGENERICTRAP)
  #define G4UGenericTrap G4GenericTrap
  #include "G4UGenericTrap.hh"
#else

#include <vector>

#include "G4TwoVector.hh"
#include "G4VSolid.hh"
#include "globals.hh"

class G4VFacet;
class G4TessellatedSolid;

class G4GenericTrap : public G4VSolid
{
  public:  // with description

     G4GenericTrap( const G4String& name, G4double halfZ,
                    const std::vector<G4TwoVector>& vertices );
       // Constructor

     ~G4GenericTrap();
       // Destructor

    // Accessors

    inline G4double    GetZHalfLength() const;
    inline G4int       GetNofVertices() const;
    inline G4TwoVector GetVertex(G4int index) const;
    inline const std::vector<G4TwoVector>& GetVertices() const;
    inline G4double    GetTwistAngle(G4int index) const;
    inline G4bool      IsTwisted() const;
    inline G4int       GetVisSubdivisions() const;
    inline void        SetVisSubdivisions(G4int subdiv);

    // Solid methods                                

    EInside Inside(const G4ThreeVector& p) const;
    G4ThreeVector SurfaceNormal(const G4ThreeVector& p) const;
    G4double DistanceToIn(const G4ThreeVector& p,
                          const G4ThreeVector& v) const;
    G4double DistanceToIn(const G4ThreeVector& p) const;
    G4double DistanceToOut(const G4ThreeVector& p,
                           const G4ThreeVector& v,
                           const G4bool calcNorm = false,
                                 G4bool *validNorm = 0,
                                 G4ThreeVector *n = 0) const;
    G4double DistanceToOut(const G4ThreeVector& p) const;
    void BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const;
    G4bool CalculateExtent(const EAxis pAxis,
                           const G4VoxelLimits& pVoxelLimit,
                           const G4AffineTransform& pTransform,
                                 G4double& pmin, G4double& pmax) const;

    G4GeometryType GetEntityType() const;

    G4VSolid* Clone() const;

    std::ostream& StreamInfo(std::ostream& os) const;

    G4ThreeVector GetPointOnSurface() const ;

    G4double GetCubicVolume();
    G4double GetSurfaceArea();

    // Visualisation functions
  
    G4Polyhedron* GetPolyhedron () const;
    void DescribeYourselfTo(G4VGraphicsScene& scene) const;
    G4VisExtent   GetExtent() const;
    G4Polyhedron* CreatePolyhedron() const;
       
  public:  // without description

    G4GenericTrap(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4GenericTrap(const G4GenericTrap& rhs);
    G4GenericTrap& operator=(const G4GenericTrap& rhs); 
      // Copy constructor and assignment operator.

  private:

    // Internal methods

    inline void SetTwistAngle(G4int index, G4double twist);
    G4bool  ComputeIsTwisted() ;
    G4bool  CheckOrder(const std::vector<G4TwoVector>& vertices) const;
    G4bool  IsSegCrossing(const G4TwoVector& a, const G4TwoVector& b, 
                          const G4TwoVector& c, const G4TwoVector& d) const;
    G4bool  IsSegCrossingZ(const G4TwoVector& a, const G4TwoVector& b, 
                           const G4TwoVector& c, const G4TwoVector& d) const;
    void ReorderVertices(std::vector<G4ThreeVector>& vertices) const;
    void ComputeBBox();
    inline G4ThreeVector GetMinimumBBox() const;
    inline G4ThreeVector GetMaximumBBox() const;
      
    G4VFacet* MakeDownFacet(const std::vector<G4ThreeVector>& fromVertices, 
                            G4int ind1, G4int ind2, G4int ind3) const;
    G4VFacet* MakeUpFacet(const std::vector<G4ThreeVector>& fromVertices, 
                            G4int ind1, G4int ind2, G4int ind3) const;      
    G4VFacet* MakeSideFacet(const G4ThreeVector& downVertex0, 
                            const G4ThreeVector& downVertex1,
                            const G4ThreeVector& upVertex1,
                            const G4ThreeVector& upVertex0) const;
    G4TessellatedSolid* CreateTessellatedSolid() const;
     
    EInside InsidePolygone(const G4ThreeVector& p,
                           const std::vector<G4TwoVector>& poly) const;
    G4double DistToPlane(const G4ThreeVector& p,
                         const G4ThreeVector& v, const G4int ipl) const ;
    G4double DistToTriangle(const G4ThreeVector& p,
                            const G4ThreeVector& v, const G4int ipl) const;
    G4ThreeVector NormalToPlane(const G4ThreeVector& p,
                                const G4int ipl) const;
    G4double SafetyToFace(const G4ThreeVector& p, const G4int iseg) const;
    G4double GetFaceSurfaceArea(const G4ThreeVector& p0,
                                const G4ThreeVector& p1,
                                const G4ThreeVector& p2,
                                const G4ThreeVector& p3) const;
    G4double GetFaceCubicVolume(const G4ThreeVector& p0,
                                const G4ThreeVector& p1,
                                const G4ThreeVector& p2,
                                const G4ThreeVector& p3) const;
  protected:

     mutable G4bool fRebuildPolyhedron;
     mutable G4Polyhedron*   fpPolyhedron;

  private:

    // static data members

    static const G4int       fgkNofVertices;
    static const G4double    fgkTolerance;

    G4double halfCarTolerance;

    // data members

    G4double                 fDz;
    std::vector<G4TwoVector> fVertices;
    G4bool                   fIsTwisted;
    G4double                 fTwist[4];
    G4TessellatedSolid*      fTessellatedSolid;
    G4ThreeVector            fMinBBoxVector;
    G4ThreeVector            fMaxBBoxVector;
    G4int                    fVisSubdivisions;

    enum ESide {kUndefined,kXY0,kXY1,kXY2,kXY3,kMZ,kPZ};
      // Codes for faces (kXY[num]=num of lateral face,kMZ= minus z face etc)

    G4double                 fSurfaceArea;
    G4double                 fCubicVolume;
      // Surface and Volume
};    

#include "G4GenericTrap.icc"

#endif

#endif
