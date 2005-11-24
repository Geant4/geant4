// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		G4TessellatedSolid.hh
//
// Date:		15/06/2005
// Author:		P R Truscott
// Organisation:	QinetiQ Ltd, UK
// Customer:		UK Ministry of Defence : RAO CRP TD Electronic Systems
// Contract:		C/MAT/N03517
//
// This software is the intelectual property of QinetiQ Ltd, subject
// DEFCON 705 IPR conditions.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
// 22 November 2005, F Lei, 
//  - Added GetPolyhedron()
//
// 31 October 2004, P R Truscott, QinetiQ Ltd, UK
// Created.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DISCLAIMER
// ----------
//
//
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DESCRIPTION
// -----------
//    G4TessellatedSolid is a special Geant4 Solid defined by a number of 
//    G4VFacet.It is important that the supplied facets shall form a fully 
//    enclose space which is the solid. 
//    At the moment only two types of facet can be used for the construction of 
//    a G4TessellatedSolid, i.e. the G4TriangularFacet and G4QuadrangularFacet.
//
//    How to contruct a G4TessellatedSolid
//  
//       .....
////      First declare a tessellated solid 
//           G4TessellatedSolid solidTarget = new G4TessellatedSolid("Solid_name");
////      Define the facets which form the solid
// 
//           G4double targetSiz = 10*cm ;
//           G4TriangularFacet *facet1 = new
//           G4TriangularFacet (G4ThreeVector(-targetSize,-targetSize,        0.0),
//                     G4ThreeVector(+targetSize,-targetSize,        0.0),
//                     G4ThreeVector(        0.0,        0.0,+targetSize),
//                     ABSOLUTE);
//           G4TriangularFacet *facet2 = new
//           G4TriangularFacet (G4ThreeVector(+targetSize,-targetSize,        0.0),
//                                G4ThreeVector(+targetSize,+targetSize,        0.0),
//                                G4ThreeVector(        0.0,        0.0,+targetSize),
//                                ABSOLUTE);
//           G4TriangularFacet *facet3 = new
//           G4TriangularFacet (G4ThreeVector(+targetSize,+targetSize,        0.0),
//                                G4ThreeVector(-targetSize,+targetSize,        0.0),
//                                G4ThreeVector(        0.0,        0.0,+targetSize),
//                                ABSOLUTE);
//           G4TriangularFacet *facet4 = new
//           G4TriangularFacet (G4ThreeVector(-targetSize,+targetSize,        0.0),
//                                G4ThreeVector(-targetSize,-targetSize,        0.0),
//                                G4ThreeVector(        0.0,        0.0,+targetSize),
//                                ABSOLUTE);
//           G4QuadrangularFacet *facet5 = new
//           G4QuadrangularFacet (G4ThreeVector(-targetSize,-targetSize,        0.0),
//                                G4ThreeVector(-targetSize,+targetSize,        0.0),
//                                G4ThreeVector(+targetSize,+targetSize,        0.0),
//                                G4ThreeVector(+targetSize,-targetSize,        0.0),
//                                ABSOLUTE);
////      Noew add the facets to the solid     
//             solidTarget->AddFacet((G4VFacet*) facet1);
//             solidTarget->AddFacet((G4VFacet*) facet2);
//             solidTarget->AddFacet((G4VFacet*) facet3);
//             solidTarget->AddFacet((G4VFacet*) facet4);
//             solidTarget->AddFacet((G4VFacet*) facet5);
////      Finally declare the solid is complete  
//             solidTarget->SetSolidClosed(true);
//
//  ...............
//
///////////////////////////////////////////////////////////////////////////////
//
//
#ifndef G4TessellatedSolid_hh
#define G4TessellatedSolid_hh 1

#include "G4VSolid.hh"
#include "G4VFacet.hh"
#include "G4VGraphicsScene.hh"
#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4VisExtent.hh"
#include "globals.hh"

#include <iostream>
#include <vector>
#include <map>

using namespace std;

class G4TessellatedSolid : public G4VSolid
{
  public:
    G4TessellatedSolid ();
    G4TessellatedSolid (const G4String &name);
    ~G4TessellatedSolid ();
    
    G4TessellatedSolid (const G4TessellatedSolid &s);
    const G4TessellatedSolid &operator= (const G4TessellatedSolid &s);
    const G4TessellatedSolid &operator+= (const G4TessellatedSolid &right);
    
    G4bool AddFacet (G4VFacet *aFacet);
    G4VFacet *GetFacet (size_t i) const;
    size_t GetNumberOfFacets () const;
    
//    G4double GetCubicVolume ();
    
//    void ComputeDimensions (G4VPVParameterisation* p, const G4int n,
//      const G4VPhysicalVolume* pRep) const;
    
    EInside Inside (const G4ThreeVector &p) const;
    G4ThreeVector SurfaceNormal (const G4ThreeVector &p) const;
    G4double DistanceToIn (const G4ThreeVector &p, const G4ThreeVector &v)
      const;
    G4double DistanceToIn (const G4ThreeVector &p) const;
    G4double DistanceToOut (const G4ThreeVector &p, const G4ThreeVector &v,
      const G4bool calcNorm=false, G4bool *validNorm=0,G4ThreeVector *n=0) const;
    G4double DistanceToOut (const G4ThreeVector &p) const;
    G4GeometryType GetEntityType () const;
    
    void SetSolidClosed (const G4bool t);
    G4bool GetSolidClosed () const;
        
    G4bool CalculateExtent(const EAxis pAxis, const G4VoxelLimits& pVoxelLimit,
      const G4AffineTransform& pTransform, G4double& pMin,
      G4double& pMax) const;

    std::ostream &StreamInfo(std::ostream &os) const;
  // Functions for visualization
 
    void          DescribeYourselfTo (G4VGraphicsScene& scene) const;
    G4VisExtent   GetExtent () const;
    G4double      GetMinXExtent () const;
    G4double      GetMaxXExtent () const;
    G4double      GetMinYExtent () const;
    G4double      GetMaxYExtent () const;
    G4double      GetMinZExtent () const;
    G4double      GetMaxZExtent () const;
    G4Polyhedron* CreatePolyhedron () const;
    G4Polyhedron* GetPolyhedron      () const;
    G4NURBS*      CreateNURBS () const;
 
  protected:  // with description
 
    void DeleteObjects ();
    void CopyObjects (const G4TessellatedSolid &s);
    G4ThreeVectorList*
    CreateRotatedVertices(const G4AffineTransform& pTransform) const;
      // Create the List of transformed vertices in the format required
      // for G4VSolid:: ClipCrossSection and ClipBetweenSections.

  private: 

    mutable G4Polyhedron* fpPolyhedron;

    std::vector<G4VFacet *>  facets;
    G4GeometryType           geometryType;
    G4double                 cubicVolume;
    G4ThreeVectorList        vertexList;
    G4double                 xMinExtent;
    G4double                 xMaxExtent;
    G4double                 yMinExtent;
    G4double                 yMaxExtent;
    G4double                 zMinExtent;
    G4double                 zMaxExtent;
    G4bool                   solidClosed;
    
    G4double                 dirTolerance;
    
};
#endif
///////////////////////////////////////////////////////////////////////////////
//


