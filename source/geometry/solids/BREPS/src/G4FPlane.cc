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
// $Id: G4FPlane.cc,v 1.17 2010-07-07 14:45:31 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4FPlane.cc
//
// ----------------------------------------------------------------------
// Corrections by S.Giani:
// - The constructor using iVec now properly stores both the internal and
//   external boundaries in the bounds vector.
// - Proper initialization of sameSense in both the constructors. 
// - Addition of third argument (sense) in the second constructor to ensure
//   consistent setting of the normal in all the client code.
// - Proper use of the tolerance in the Intersect function.
// ----------------------------------------------------------------------

#include "G4FPlane.hh"
#include "G4CompositeCurve.hh"


G4FPlane::G4FPlane( const G4Vector3D& direction,
		    const G4Vector3D& axis     , 
		    const G4Point3D&  Pt0, G4int sense )
  : pplace(direction, axis, Pt0), Convex(0), projectedBoundary(0)
{
  G4Point3D Pt1 = G4Point3D( Pt0 + direction );

  // The plane include direction and axis is the normal,
  // so axis^direction is included in the plane
  G4Point3D Pt2 = G4Point3D( Pt0 + axis.cross(direction) );

  G4Ray::CalcPlane3Pts( Pl, Pt0, Pt1, Pt2 );

  active   = 1;
  sameSense = sense;
  CalcNormal();
  distance = kInfinity;
  Type     = 1;
}


G4FPlane::G4FPlane(const G4Point3DVector* pVec,
                   const G4Point3DVector* iVec,
		   G4int sense)
  : pplace( (*pVec)[0]-(*pVec)[1],                    // direction
	    ((*pVec)[pVec->size()-1]-(*pVec)[0])
	    .cross((*pVec)[0]-(*pVec)[1]),            // axis
	    (*pVec)[0]                             )  // location

{
  G4Ray::CalcPlane3Pts( Pl, (*pVec)[0], (*pVec)[1], (*pVec)[2] );
 
  G4CurveVector bounds;
  G4CompositeCurve* polygon;

  projectedBoundary = new G4SurfaceBoundary;

  sameSense = sense;

  // Outer boundary

  polygon= new G4CompositeCurve(*pVec);
 
  for (size_t i=0; i< polygon->GetSegments().size(); i++) 
    polygon->GetSegments()[i]->SetSameSense(sameSense);

  bounds.push_back(polygon);
  
  // Eventual inner boundary
  
  if (iVec) 
  {
    polygon= new G4CompositeCurve(*iVec);

    for (size_t i=0; i< polygon->GetSegments().size(); i++) 
    polygon->GetSegments()[i]->SetSameSense(sameSense);

    bounds.push_back(polygon);
  }
  
  // Set sense for boundaries  
  
  for (size_t j=0; j< bounds.size(); j++) 
    bounds[j]->SetSameSense(sameSense);
  

  SetBoundaries(&bounds);
      
  CalcNormal();
  IsConvex();
  distance = kInfinity;
  Type=1;
}


G4FPlane::~G4FPlane()
{
  delete NormalX;
}


void G4FPlane::CalcBBox()
{
  // This is needed since the bounds are used for the Solid
  // bbox calculation. The bbox test is NOT performed for
  // planar surfaces.

  // Finds the bounds of the G4Plane surface iow
  // calculates the bounds for a bounding box
  // to the surface. The bounding box is used
  // for a preliminary check of intersection.
  
  bbox= new G4BoundingBox3D(surfaceBoundary.BBox().GetBoxMin(), 
			    surfaceBoundary.BBox().GetBoxMax());

}


void G4FPlane::CalcNormal()
{
/* 
  // Calc Normal for surface which is used for the projection
  // Make planes
  G4Vector3D norm;

  G4Vector3D RefDirection = pplace.GetRefDirection();
  G4Vector3D Axis = pplace.GetAxis();

  // L. Broglia : before in G4Placement
  if( RefDirection == Axis )
    norm = RefDirection;
  else
  {
    // L. Broglia : error on setY, and it`s better to use cross function
    // norm.setX( RefDirection.y() * Axis.z() -  RefDirection.z() * Axis.y() );
    // norm.setY( RefDirection.x() * Axis.z() -  RefDirection.z() * Axis.x() );
    // norm.setZ( RefDirection.x() * Axis.y() -  RefDirection.y() * Axis.x() );
       
    norm = RefDirection.cross(Axis);
  }
  
  //  const G4Point3D& tmp = pplace.GetSrfPoint();
  const G4Point3D tmp = pplace.GetLocation();
*/  

  // L. Broglia
  // The direction of the normal is the axis of his location
  // Its sense depend on the orientation of the bounded curve
  const G4Point3D tmp = pplace.GetLocation();
  G4Vector3D norm;
  G4int sense = GetSameSense();
  
  if (sense)
    norm = pplace.GetAxis();
  else
    norm = - pplace.GetAxis();

  NormalX =  new G4Ray(tmp, norm);
  NormalX->RayCheck();
  NormalX->CreatePlanes();
}


void G4FPlane::Project()
{
    // Project
    // const G4Plane& Plane1 = NormalX->GetPlane(1);
    // const G4Plane& Plane2 = NormalX->GetPlane(2);

    // probably not necessary
    // projections of the boundary should be handled by the intersection
    //    OuterBoundary->ProjectBoundaryTo2D(Plane1, Plane2, 0);
}


G4int G4FPlane::IsConvex() const
{
  return -1;  
}


G4int G4FPlane::Intersect(const G4Ray& rayref)
{
  // This function count the number of intersections of a 
  // bounded surface by a ray.
  

  // Find the intersection with the infinite plane
  Intersected =1;

  // s is solution, line is p + tq, n is G4Plane Normal, r is point on G4Plane 
  // all parameters are pointers to arrays of three elements
 
  hitpoint = PINFINITY;
  register G4double a, b, t;

  register const G4Vector3D& RayDir   = rayref.GetDir();
  register const G4Point3D&  RayStart = rayref.GetStart();

  G4double dirx =  RayDir.x();
  G4double diry =  RayDir.y();
  G4double dirz =  RayDir.z();

  G4Vector3D norm = (*NormalX).GetDir();
  G4Point3D  srf_point = pplace.GetLocation();

  b = norm.x() * dirx + norm.y() * diry + norm.z() * dirz;

  if ( std::fabs(b) < perMillion )    
  {
    // G4cout << "\nLine is parallel to G4Plane.No Hit.";
  }  
  else
  {
    G4double startx =  RayStart.x();
    G4double starty =  RayStart.y();
    G4double startz =  RayStart.z();    
    
    a = norm.x() * (srf_point.x() - startx) + 
        norm.y() * (srf_point.y() - starty) + 
        norm.z() * (srf_point.z() - startz)   ;
    
    t = a/b;
    
    // substitute t into line equation
    // to calculate final solution     
    G4double solx,soly,solz;
    solx = startx + t * dirx;
    soly = starty + t * diry;
    solz = startz + t * dirz;

    // solve tolerance problem
    if( (t*dirx >= -kCarTolerance/2) && (t*dirx <= kCarTolerance/2) )
      solx = startx;

    if( (t*diry >= -kCarTolerance/2) && (t*diry <= kCarTolerance/2) )
      soly = starty;

    if( (t*dirz >= -kCarTolerance/2) && (t*dirz <= kCarTolerance/2) )
      solz = startz;
    
    G4bool xhit = (dirx < 0 && solx <= startx) || (dirx >= 0 && solx >= startx);
    G4bool yhit = (diry < 0 && soly <= starty) || (diry >= 0 && soly >= starty);
    G4bool zhit = (dirz < 0 && solz <= startz) || (dirz >= 0 && solz >= startz);
    
    if( xhit && yhit && zhit ) {
      hitpoint= G4Point3D(solx, soly, solz);
    }
  }
   
  // closest_hit is a public Point3D in G4Surface
  closest_hit = hitpoint;
  
  if(closest_hit.x() == kInfinity)
  {
    // no hit
    active=0;
    SetDistance(kInfinity);
    return 0;
  }
  else
  {
    // calculate the squared distance from the point to the intersection
    // and set it in the distance data member (all clients know they have
    // to take the sqrt)
    SetDistance( RayStart.distance2(closest_hit) );   

    // now, we have to verify that the hit point founded
    // is included into the G4FPlane boundaries

    // project the hit to the xy plane,
    // with the same projection that took the boundary
    // into projectedBoundary
    G4Point3D projectedHit= pplace.GetToPlacementCoordinates() * closest_hit;
    
    // test ray from the hit on the xy plane
    G4Ray testRay( projectedHit, G4Vector3D(1, 0.01, 0) );

    // check if it intersects the boundary
    G4int nbinter = projectedBoundary->IntersectRay2D(testRay);

    // If this number is par, it`s signify that the projected point  
    // is outside the projected surface, so the hit point is outside
    // the bounded surface
    if(nbinter&1)
    {
      // the intersection point is into the boundaries
      // check if the intersection point is on the surface
      if(distance <= kCarTolerance*0.5*kCarTolerance*0.5)
      {
	// the point is on the surface, set the distance to 0            
	SetDistance(0);         
      }
      else
      {
	// the point is outside the surface
      }
      
      return 1 ;      
    }
    else
    {
      // the intersection point is out the boundaries
      // it is not a real intersection
      active=0;
      SetDistance(kInfinity);
      return 0;
    }
  }
}


G4double G4FPlane::ClosestDistanceToPoint(const G4Point3D& Pt)
{
  // Calculates signed distance of point Pt to G4Plane Pl
  // Be careful, the equation of the plane is :
  // ax + by + cz = d
  G4double dist = Pt.x()*Pl.a + Pt.y()*Pl.b + Pt.z()*Pl.c - Pl.d;

  return dist;
}


void G4FPlane::InitBounded()
{
  // L. Broglia

  projectedBoundary =
    surfaceBoundary.Project( pplace.GetToPlacementCoordinates() );
}

G4double G4FPlane::HowNear( const G4Vector3D& Pt ) const
{
  G4double hownear = Pt.x()*Pl.a + Pt.y()*Pl.b + Pt.z()*Pl.c - Pl.d;

  return hownear;
}
