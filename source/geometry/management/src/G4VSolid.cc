// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VSolid.cc,v 1.5 2000-11-01 15:39:36 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4VSolid
//
// Implementation for solid base class
//
//
// History:
//  10.07.95 P.Kent Added == operator, solid Store entry
//  30.06.95 P.Kent

#include "G4VSolid.hh"
#include "G4SolidStore.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4VisExtent.hh"

// Constructor
//  - Copies name
//  - Add ourselves to solid Store
G4VSolid::G4VSolid(const G4String& name) :
fshapeName(name) 
{
    G4SolidStore::GetInstance()->append(this);
}

// Destructor (virtual)
// - Remove ourselves from solid Store
G4VSolid::~G4VSolid()
{
    G4SolidStore::GetInstance()->remove(this);
}
	
// Returns name by value
G4String G4VSolid::GetName() const
{
    return fshapeName;
}	

void G4VSolid::SetName(const G4String& name)
{
    fshapeName=name;
}	

// Throw exception if ComputeDimensions called for illegal derived class
void G4VSolid::ComputeDimensions(G4VPVParameterisation* p,
	                         const G4int n,
                                 const G4VPhysicalVolume* pRep)
{
    G4Exception("G4VSolid::ComputeDimensions called illegally: not overloaded by derived class");
}

// Calculate the maximum and minimum extents of the convex polygon pPolygon
// along the axis pAxis, within the limits pVoxelLimit
void G4VSolid::CalculateClippedPolygonExtent(G4ThreeVectorList& pPolygon,
					  const G4VoxelLimits& pVoxelLimit,
					  const EAxis pAxis, 
					  G4double& pMin, G4double& pMax) const
{
    G4int noLeft,i;
    G4double component;
    ClipPolygon(pPolygon,pVoxelLimit);
    noLeft=pPolygon.entries();
    if (noLeft)
	{
	    for (i=0;i<noLeft;i++)
		{
		    component=pPolygon(i).operator()(pAxis);
		    if (component<pMin)
			{
			    pMin=component;
			}
		    else if (component>pMax)
			{
			    pMax=component;
			}
		}
	}
}
 
// Calculate the maximum and minimum extents of the polygon described
// by the vertices: pSectionIndex->pSectionIndex+1->
//                   pSectionIndex+2->pSectionIndex+3->pSectionIndex
// in the List pVertices
//
// If the minimum is <pMin pMin is set to the new minimum
// If the maximum is >pMax pMax is set to the new maximum
//
// No modifications are made to pVertices
void G4VSolid::ClipCrossSection(G4ThreeVectorList* pVertices,
			     const G4int pSectionIndex,
			     const G4VoxelLimits& pVoxelLimit,
			     const EAxis pAxis, 
			     G4double& pMin, G4double& pMax) const
{

    G4ThreeVectorList polygon;
    polygon.append(pVertices->operator()(pSectionIndex));
    polygon.append(pVertices->operator()(pSectionIndex+1));
    polygon.append(pVertices->operator()(pSectionIndex+2));
    polygon.append(pVertices->operator()(pSectionIndex+3));
    CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
    return;
}

// Calculate the maximum and minimum extents of the polygons
// joining the CrossSections at pSectionIndex->pSectionIndex+3 and
//                              pSectionIndex+4->pSectionIndex7
//
// in the List pVertices, within the boundaries of the voxel limits pVoxelLimit
//
// If the minimum is <pMin pMin is set to the new minimum
// If the maximum is >pMax pMax is set to the new maximum
//
// No modifications are made to pVertices

void G4VSolid::ClipBetweenSections(G4ThreeVectorList* pVertices,
			     const G4int pSectionIndex,
			     const G4VoxelLimits& pVoxelLimit,
			     const EAxis pAxis, 
			     G4double& pMin, G4double& pMax) const
{
    G4ThreeVectorList polygon;
    polygon.append(pVertices->operator()(pSectionIndex));
    polygon.append(pVertices->operator()(pSectionIndex+4));
    polygon.append(pVertices->operator()(pSectionIndex+5));
    polygon.append(pVertices->operator()(pSectionIndex+1));
    CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
    polygon.clear();
    polygon.append(pVertices->operator()(pSectionIndex+1));
    polygon.append(pVertices->operator()(pSectionIndex+5));
    polygon.append(pVertices->operator()(pSectionIndex+6));
    polygon.append(pVertices->operator()(pSectionIndex+2));
    CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
    polygon.clear();
    polygon.append(pVertices->operator()(pSectionIndex+2));
    polygon.append(pVertices->operator()(pSectionIndex+6));
    polygon.append(pVertices->operator()(pSectionIndex+7));
    polygon.append(pVertices->operator()(pSectionIndex+3));
    CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
    polygon.clear();
    polygon.append(pVertices->operator()(pSectionIndex+3));
    polygon.append(pVertices->operator()(pSectionIndex+7));
    polygon.append(pVertices->operator()(pSectionIndex+4));
    polygon.append(pVertices->operator()(pSectionIndex));
    CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
    return;
}

// Clip the convex polygon described by the vertices at
// pSectionIndex ->pSectionIndex+3 within pVertices to the limits pVoxelLimit
//
// Set pMin to the smallest
//
// Calculate the extent of the polygon along pAxis, when clipped to the
// limits pVoxelLimit. If the polygon exists after clippin, set pMin to
// the polygon's minimum extent along the axis if <pMin, and set pMax to
// the polygon's maximum extent along the axis if >pMax.
//
// The polygon is described by a set of vectors, where each vector represents
// a vertex, so that the polygon is described by the vertex sequence:
//   0th->1st 1st->2nd 2nd->... nth->0th
//
// Modifications to the polygon are made
//
// NOTE: Execessive copying during clipping
void G4VSolid::ClipPolygon(G4ThreeVectorList& pPolygon,
			const G4VoxelLimits& pVoxelLimit) const
{
    G4ThreeVectorList outputPolygon;
    if (pVoxelLimit.IsLimited())
	{
	    if (pVoxelLimit.IsXLimited())
		{
		    G4VoxelLimits simpleLimit1;
		    simpleLimit1.AddLimit(kXAxis,pVoxelLimit.GetMinXExtent(),kInfinity);
		    ClipPolygonToSimpleLimits(pPolygon,outputPolygon,
					      simpleLimit1);
		    pPolygon.clear();
		    if (!outputPolygon.entries())
			{
			    return;
			}

		    G4VoxelLimits simpleLimit2;
		    simpleLimit2.AddLimit(kXAxis,-kInfinity,pVoxelLimit.GetMaxXExtent());
		    ClipPolygonToSimpleLimits(outputPolygon,pPolygon,simpleLimit2);
		    if (!pPolygon.entries())
			{
			    return;
			}
		    else
			{
			    outputPolygon.clear();
			}
		}

	    if (pVoxelLimit.IsYLimited())
		{
		    G4VoxelLimits simpleLimit1;
		    simpleLimit1.AddLimit(kYAxis,pVoxelLimit.GetMinYExtent(),kInfinity);
		    ClipPolygonToSimpleLimits(pPolygon,outputPolygon,
					      simpleLimit1);
// Must always clear pPolygon - for clip to simpleLimit2 and incase of
// early exit
		    pPolygon.clear();
		    if (!outputPolygon.entries())
			{
			    return;
			}

		    G4VoxelLimits simpleLimit2;
		    simpleLimit2.AddLimit(kYAxis,-kInfinity,pVoxelLimit.GetMaxYExtent());
		    ClipPolygonToSimpleLimits(outputPolygon,pPolygon,simpleLimit2);
		    if (!pPolygon.entries())
			{
			    return;
			}
		    else
			{
			    outputPolygon.clear();
			}
		}

	    if (pVoxelLimit.IsZLimited())
		{
		    G4VoxelLimits simpleLimit1;
		    simpleLimit1.AddLimit(kZAxis,pVoxelLimit.GetMinZExtent(),kInfinity);
		    ClipPolygonToSimpleLimits(pPolygon,outputPolygon,
					      simpleLimit1);
// Must always clear pPolygon - for clip to simpleLimit2 and incase of
// early exit
		    pPolygon.clear();
		    if (!outputPolygon.entries())
			{
			    return;
			}

		    G4VoxelLimits simpleLimit2;
		    simpleLimit2.AddLimit(kZAxis,-kInfinity,pVoxelLimit.GetMaxZExtent());
		    ClipPolygonToSimpleLimits(outputPolygon,pPolygon,simpleLimit2);

// Return after final clip - no cleanup
		}
	}
}

// pVoxelLimits must be only limited along one axis, and either the maximum
// along the axis must be +kInfinity, or the minimum -kInfinity
void G4VSolid::ClipPolygonToSimpleLimits(G4ThreeVectorList& pPolygon,
				      G4ThreeVectorList& outputPolygon,
				      const G4VoxelLimits& pVoxelLimit) const
{
    G4int i;
    G4int noVertices=pPolygon.entries();
    G4ThreeVector vEnd,vStart;
    for (i=0;i<noVertices;i++)
	{
	    vStart=pPolygon(i);
	    if (i==noVertices-1)
		{
		    vEnd=pPolygon(0);
		}
	    else
		{
		    vEnd=pPolygon(i+1);
		}

	    if (pVoxelLimit.Inside(vStart))
		{
		    if (pVoxelLimit.Inside(vEnd))
			{
// vStart and vEnd inside -> output end point
			    outputPolygon.insert(vEnd);
			}
		    else
			{
// vStart inside, vEnd outside -> output crossing point
			    pVoxelLimit.ClipToLimits(vStart,vEnd);
			    outputPolygon.insert(vEnd);
			}
		    
		}
	    else
		{
		    if (pVoxelLimit.Inside(vEnd))
			{
// vStart outside, vEnd inside -> output inside section
			    pVoxelLimit.ClipToLimits(vStart,vEnd);
			    outputPolygon.insert(vStart);
			    outputPolygon.insert(vEnd);
			}
		    else
// Both point outside -> no output
			{
			}

		}
	}

}

 const G4VSolid* G4VSolid::GetConstituentSolid(G4int no) const
{ return 0; } 

       G4VSolid* G4VSolid::GetConstituentSolid(G4int no)
{ return 0; } 

 const G4DisplacedSolid* G4VSolid::GetDisplacedSolidPtr() const
{ return 0; } 

       G4DisplacedSolid* G4VSolid::GetDisplacedSolidPtr() 
{ return 0; } 

G4VisExtent G4VSolid::GetExtent () const {
  G4VisExtent extent;
  G4VoxelLimits voxelLimits;  // Defaults to "infinite" limits.
  G4AffineTransform affineTransform;
  G4double vmin, vmax;
  CalculateExtent(kXAxis,voxelLimits,affineTransform,vmin,vmax);
  extent.SetXmin (vmin);
  extent.SetXmax (vmax);
  CalculateExtent(kYAxis,voxelLimits,affineTransform,vmin,vmax);
  extent.SetYmin (vmin);
  extent.SetYmax (vmax);
  CalculateExtent(kZAxis,voxelLimits,affineTransform,vmin,vmax);
  extent.SetZmin (vmin);
  extent.SetZmax (vmax);
  return extent;
}

G4Polyhedron* G4VSolid::CreatePolyhedron () const
{
  return 0;
}

G4NURBS* G4VSolid::CreateNURBS () const
{
  return 0;
}
