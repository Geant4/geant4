//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//

//
// class G4VSolid
//
// Implementation for solid base class
//
// History:
//
//  06.12.02 V.Grichine, bug fixed in ClipPoligon: clip as P.Kent version
//                       
//  10.05.02 V.Grichine, bug fixed in ClipPoligon: clip only other axis and limited
//                       voxels
//  15.04.02 V.Grichine, bug fixed in ClipPoligon: clip only one axis
//  13.03.02 V.Grichine, cosmetics of voxel limit functions  
//  15.11.00 D.Williams, V.Grichine change in CalculateClippedPolygonExtent:
//                                  else if(component>pMax) ---> if
//  10.07.95 P.Kent Added == operator, solid Store entry
//  30.06.95 P.Kent
//                                  

#include "G4VSolid.hh"
#include "G4SolidStore.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4VisExtent.hh"

//////////////////////////////////////////////////////////////////////////
//
// Constructor
//  - Copies name
//  - Add ourselves to solid Store
G4VSolid::G4VSolid(const G4String& name)
  : fshapeName(name) 
{
    G4SolidStore::GetInstance()->Register(this);
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor (virtual)
// - Remove ourselves from solid Store
G4VSolid::~G4VSolid()
{
    G4SolidStore::GetInstance()->DeRegister(this);
}

//////////////////////////////////////////////////////////////////////////
//
// Streaming operator dumping solid contents
G4std::ostream& operator<< ( G4std::ostream& os, const G4VSolid& e )
{
    return e.StreamInfo(os);
}

//////////////////////////////////////////////////////////////////////////
//
// Throw exception if ComputeDimensions called for illegal derived class

void G4VSolid::ComputeDimensions(G4VPVParameterisation* p,
	                         const G4int n,
                                 const G4VPhysicalVolume* pRep)
{
    G4Exception("G4VSolid::ComputeDimensions called illegally: not overloaded by derived class");
}

///////////////////////////////////////////////////////////////////////////
// 
// Calculate the maximum and minimum extents of the polygon described
// by the vertices: pSectionIndex->pSectionIndex+1->
//                   pSectionIndex+2->pSectionIndex+3->pSectionIndex
// in the List pVertices
//
// If the minimum is <pMin pMin is set to the new minimum
// If the maximum is >pMax pMax is set to the new maximum
//
// No modifications are made to pVertices
//

void G4VSolid::ClipCrossSection(       G4ThreeVectorList* pVertices,
			         const G4int pSectionIndex,
			         const G4VoxelLimits& pVoxelLimit,
			         const EAxis pAxis, 
			               G4double& pMin, G4double& pMax) const
{

  G4ThreeVectorList polygon;
  polygon.push_back((*pVertices)[pSectionIndex]);
  polygon.push_back((*pVertices)[pSectionIndex+1]);
  polygon.push_back((*pVertices)[pSectionIndex+2]);
  polygon.push_back((*pVertices)[pSectionIndex+3]);
  //  G4cout<<"ClipCrossSection: 0-1-2-3"<<G4endl;
  CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
  return;
}

//////////////////////////////////////////////////////////////////////////////////
//
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

void G4VSolid::ClipBetweenSections(      G4ThreeVectorList* pVertices,
			           const G4int pSectionIndex,
			           const G4VoxelLimits& pVoxelLimit,
			           const EAxis pAxis, 
			                 G4double& pMin, G4double& pMax) const
{
  G4ThreeVectorList polygon;
  polygon.push_back((*pVertices)[pSectionIndex]);
  polygon.push_back((*pVertices)[pSectionIndex+4]);
  polygon.push_back((*pVertices)[pSectionIndex+5]);
  polygon.push_back((*pVertices)[pSectionIndex+1]);
  // G4cout<<"ClipBetweenSections: 0-4-5-1"<<G4endl;
  CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
  polygon.clear();

  polygon.push_back((*pVertices)[pSectionIndex+1]);
  polygon.push_back((*pVertices)[pSectionIndex+5]);
  polygon.push_back((*pVertices)[pSectionIndex+6]);
  polygon.push_back((*pVertices)[pSectionIndex+2]);
  // G4cout<<"ClipBetweenSections: 1-5-6-2"<<G4endl;
  CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
  polygon.clear();

  polygon.push_back((*pVertices)[pSectionIndex+2]);
  polygon.push_back((*pVertices)[pSectionIndex+6]);
  polygon.push_back((*pVertices)[pSectionIndex+7]);
  polygon.push_back((*pVertices)[pSectionIndex+3]);
  //  G4cout<<"ClipBetweenSections: 2-6-7-3"<<G4endl;
  CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
  polygon.clear();

  polygon.push_back((*pVertices)[pSectionIndex+3]);
  polygon.push_back((*pVertices)[pSectionIndex+7]);
  polygon.push_back((*pVertices)[pSectionIndex+4]);
  polygon.push_back((*pVertices)[pSectionIndex]);
  //  G4cout<<"ClipBetweenSections: 3-7-4-0"<<G4endl;
  CalculateClippedPolygonExtent(polygon,pVoxelLimit,pAxis,pMin,pMax);
  return;
}


///////////////////////////////////////////////////////////////////////////////
//
// Calculate the maximum and minimum extents of the convex polygon pPolygon
// along the axis pAxis, within the limits pVoxelLimit
//

void G4VSolid::CalculateClippedPolygonExtent(G4ThreeVectorList& pPolygon,
					  const G4VoxelLimits& pVoxelLimit,
					  const EAxis pAxis, 
					  G4double& pMin, G4double& pMax) const
{
  G4int noLeft,i;
  G4double component;
  /*  
  G4cout<<G4endl;
  for(i = 0 ; i < pPolygon.size() ; i++ )
  {
      G4cout<<i<<"\t"<<"p.x = "<<pPolygon[i].operator()(pAxis)<<"\t"
	//            <<"p.y = "<<pPolygon[i].y()<<"\t"
        //             <<"p.z = "<<pPolygon[i].z()<<"\t"
            <<G4endl;
  }	    
  G4cout<<G4endl;
  */  
  ClipPolygon(pPolygon,pVoxelLimit,pAxis);
  noLeft = pPolygon.size();

  if ( noLeft )
  {
    //  G4cout<<G4endl;
    for (i=0;i<noLeft;i++)
    {
      component = pPolygon[i].operator()(pAxis);
      //  G4cout <<i<<"\t"<<component<<G4endl;
 
      if (component < pMin) 
      { 
	//  G4cout <<i<<"\t"<<"Pmin = "<<component<<G4endl;
        pMin = component;      
      }
      if (component > pMax)
      {  
	//  G4cout <<i<<"\t"<<"PMax = "<<component<<G4endl;
        pMax = component;  
      }    
    }
    //  G4cout<<G4endl;
  }
  // G4cout<<"pMin = "<<pMin<<"\t"<<"pMax = "<<pMax<<G4endl;
}

/////////////////////////////////////////////////////////////////////////////
//
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

void G4VSolid::ClipPolygon(      G4ThreeVectorList& pPolygon,
			   const G4VoxelLimits& pVoxelLimit,
                           const EAxis              pAxis          ) const
{
  G4ThreeVectorList outputPolygon;

  if ( pVoxelLimit.IsLimited() )
  {
    if (pVoxelLimit.IsXLimited() ) // && pAxis != kXAxis)
    {
      G4VoxelLimits simpleLimit1;
      simpleLimit1.AddLimit(kXAxis,pVoxelLimit.GetMinXExtent(),kInfinity);
      //  G4cout<<"MinXExtent()"<<G4endl;
      ClipPolygonToSimpleLimits(pPolygon,outputPolygon,simpleLimit1);
   
      pPolygon.clear();

      if ( !outputPolygon.size() )  return;
			
      G4VoxelLimits simpleLimit2;
      //  G4cout<<"MaxXExtent()"<<G4endl;
      simpleLimit2.AddLimit(kXAxis,-kInfinity,pVoxelLimit.GetMaxXExtent());
      ClipPolygonToSimpleLimits(outputPolygon,pPolygon,simpleLimit2);

      if ( !pPolygon.size() )       return;		
      else                          outputPolygon.clear();	
    }
	if ( pVoxelLimit.IsYLimited() ) // && pAxis != kYAxis)
    {
      G4VoxelLimits simpleLimit1;
      simpleLimit1.AddLimit(kYAxis,pVoxelLimit.GetMinYExtent(),kInfinity);
      ClipPolygonToSimpleLimits(pPolygon,outputPolygon,simpleLimit1);

// Must always clear pPolygon - for clip to simpleLimit2 and in case of
// early exit

      pPolygon.clear();

      if ( !outputPolygon.size() )  return;
			
      G4VoxelLimits simpleLimit2;
      simpleLimit2.AddLimit(kYAxis,-kInfinity,pVoxelLimit.GetMaxYExtent());
      ClipPolygonToSimpleLimits(outputPolygon,pPolygon,simpleLimit2);

      if ( !pPolygon.size() )       return;
      else                          outputPolygon.clear();	
    }
    if ( pVoxelLimit.IsZLimited() ) // && pAxis != kZAxis)
    {
      G4VoxelLimits simpleLimit1;
      simpleLimit1.AddLimit(kZAxis,pVoxelLimit.GetMinZExtent(),kInfinity);
      ClipPolygonToSimpleLimits(pPolygon,outputPolygon,simpleLimit1);

// Must always clear pPolygon - for clip to simpleLimit2 and in case of
// early exit

      pPolygon.clear();

      if ( !outputPolygon.size() )  return;
			
      G4VoxelLimits simpleLimit2;
      simpleLimit2.AddLimit(kZAxis,-kInfinity,pVoxelLimit.GetMaxZExtent());
      ClipPolygonToSimpleLimits(outputPolygon,pPolygon,simpleLimit2);

// Return after final clip - no cleanup
    }
  }
}

////////////////////////////////////////////////////////////////////////////
//
// pVoxelLimits must be only limited along one axis, and either the maximum
// along the axis must be +kInfinity, or the minimum -kInfinity

void G4VSolid::ClipPolygonToSimpleLimits( G4ThreeVectorList& pPolygon,
				          G4ThreeVectorList& outputPolygon,
				    const G4VoxelLimits& pVoxelLimit       ) const
{
  G4int i;
  G4int noVertices=pPolygon.size();
  G4ThreeVector vEnd,vStart;

  for (i = 0 ; i < noVertices ; i++ )
  {
    vStart = pPolygon[i];
    //  G4cout<<"i = "<<i<<G4endl;
    if ( i == noVertices-1 )    vEnd = pPolygon[0];
    else                        vEnd = pPolygon[i+1];
		
    if ( pVoxelLimit.Inside(vStart) )
    {
      if (pVoxelLimit.Inside(vEnd))
      {
// vStart and vEnd inside -> output end point

        outputPolygon.push_back(vEnd);
      }
      else
      {
// vStart inside, vEnd outside -> output crossing point
	// G4cout<<"vStart inside, vEnd outside"<<G4endl;
        pVoxelLimit.ClipToLimits(vStart,vEnd);
	outputPolygon.push_back(vEnd);
      }	    
    }
    else
    {
      if (pVoxelLimit.Inside(vEnd))
      {
// vStart outside, vEnd inside -> output inside section
	//  G4cout<<"vStart outside, vEnd inside"<<G4endl;

        pVoxelLimit.ClipToLimits(vStart,vEnd);
	outputPolygon.push_back(vStart);
	outputPolygon.push_back(vEnd);  
      }
      else  // Both point outside -> no output
      {
	// outputPolygon.push_back(vStart);
	// outputPolygon.push_back(vEnd);  
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

G4VisExtent G4VSolid::GetExtent () const 
{
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
