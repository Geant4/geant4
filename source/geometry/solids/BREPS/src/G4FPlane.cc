// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FPlane.cc,v 1.2 1999-01-14 16:11:19 broglia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//


#include "G4FPlane.hh"
#include "G4CompositeCurve.hh"


G4FPlane::G4FPlane( const G4Vector3D& direction,
		    const G4Vector3D& axis     , 
		    const G4Point3D&  Pt0        ):pplace(direction, axis, Pt0)
{
  G4Point3D Pt1 = Pt0 + direction;

  // The plane include direction and axis is the normal,
  // so axis^direction is included in the plane
  G4Point3D Pt2 = Pt0 + axis.cross(direction);

  G4Ray::CalcPlane3Pts( Pl, Pt0, Pt1, Pt2 );

  active   = 1;
  CalcNormal();
  distance = kInfinity;
  Type     = 1;
}


G4FPlane::G4FPlane(const G4Point3DVector* pVec, const G4Point3DVector* iVec)
  : pplace( (*pVec)[0]-(*pVec)[1],                    // direction
	    ((*pVec)[pVec->length()-1]-(*pVec)[0])
	    .cross((*pVec)[0]-(*pVec)[1]),            // axis
	    (*pVec)[0]                             )  // location

{
  G4Ray::CalcPlane3Pts( Pl, (*pVec)[0], (*pVec)[1], (*pVec)[2] );
 
  G4CurveVector bounds;
  G4CompositeCurve* polygon;

  projectedBoundary = new G4SurfaceBoundary;

  // Calcul sense of the plane, for G4BREPSolidPolyhedra 
  // rule : if direction.z > 0, sense = 1
  if(pplace.GetRefDirection().z() < 0)
    sameSense = 1; //0 
  else 
    sameSense = 1;

  polygon= new G4CompositeCurve(*pVec);

  if (iVec) 
  {
    polygon= new G4CompositeCurve(*iVec);
    bounds.insert(polygon);
  }
 
  for (G4int i=0; i< polygon->GetSegments().length(); i++) 
    polygon->GetSegments()[i]->SetSameSense(sameSense);

  bounds.insert(polygon);
  for (G4int j=0; j< bounds.length(); j++) 
    bounds[j]->SetSameSense(sameSense);
  

  SetBoundaries(&bounds);
      
  CalcNormal();
  IsConvex();
  distance = kInfinity;
  Type=1;
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
    const G4Plane& Plane1 = NormalX->GetPlane(1);
    const G4Plane& Plane2 = NormalX->GetPlane(2);

    // probably not necessary
    // projections of the boundary should be handled by the intersection
    //    OuterBoundary->ProjectBoundaryTo2D(Plane1, Plane2, 0);
}


int G4FPlane::IsConvex()
{
  return -1;  
}


int G4FPlane::Intersect(const G4Ray& rayref)
{
  Intersected =1;

  // closest_hit = pplace.EvaluateIntersection(rayref);
  // L. Broglia : before in G4Placement

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

  if ( fabs(b) < 0.001 )    
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

    if( (solx > -kCarTolerance/2) && (solx < kCarTolerance/2) )
      solx = 0;

    if( (soly > -kCarTolerance/2) && (soly < kCarTolerance/2) )
      soly = 0;

    if( (solz > -kCarTolerance/2) && (solz < kCarTolerance/2) )
      solz = 0;
    
    if(((dirx < 0 && solx < startx)||(dirx >= 0 && solx >= startx))&&
       ((diry < 0 && soly < starty)||(diry >= 0 && soly >= starty))&&
       ((dirz < 0 && solz < startz)||(dirz >= 0 && solz >= startz)))
      hitpoint= G4Point3D(solx,soly, solz);

  }
   
  // closest_hit is a public Point3D in G4Surface
  closest_hit = hitpoint;
  
  if(closest_hit.x() == kInfinity)
  {
    active=0;
    Distance(kInfinity);
    return 0;
  }
  else
  {
    Distance( RayStart.distance2(closest_hit) );

    G4Point3D hit = closest_hit;    

    // project the hit to the xy plane,
    // with the same projection that took the boundary
    // into projectedBoundary
    G4Point3D projectedHit= pplace.GetToPlacementCoordinates() * hit;
    
    // test ray from the hit on the xy plane
    // check if it intersects the boundary
    G4Ray testRay(projectedHit, G4Vector3D(1, 0.01, 0));

    G4int nbinter = projectedBoundary->IntersectRay2D(testRay);

    if(nbinter&1)
    {
      // the intersection point is into the boundaries
      // check if the intersection point is on the surface
      if(distance < kCarTolerance*0.5)
      {
	// the point is on the surface            
	Distance(0);         
      }

      return 1 ;      
    }
    else
    {
      // the intersection point is out the boundaries
      active=0;
      Distance(kInfinity);
      return 0;
    }

    
   

/*
    // if not, we are outside
    if ( is.GetDistance() >= kInfinity ) 
    {
      active=0;
      Distance(kInfinity);
      return 0;
    }


   
    // if yes, we have to check on which side of the intersected 
    // curve the hit lies
    G4Vector3D tangent;
    
    projectedBoundary->Tangent(is, tangent);
    
    // L. Broglia
    // Now replace tangent into the pplace
    tangent = pplace.GetFromPlacementCoordinates() * tangent;

    // (let's assume that the tangent is defined)
    // criterion for outside: (d x t).z() < 0
    // d = hit - is  &  t = tangent
    
    G4Point3D  Is =  pplace.GetFromPlacementCoordinates() * (is.GetPoint());
    G4Vector3D d  = hit - Is;
   
    if ( (d.cross(tangent)).z() < 0 )          
    {
      active=0;
      Distance(kInfinity);
      return 0;
    }
    
    // a real intersection point
    return 1;
*/
  }
}


G4double G4FPlane::ClosestDistanceToPoint(const G4Point3D& Pt)
{
  // Calculates signed distance of point Pt to G4Plane Pl
  // Be careful, the equation of the plane is :
  // ax + by + cz = d
  return ( Pt.x()*Pl.a + Pt.y()*Pl.b + Pt.z()*Pl.c - Pl.d);
}


void G4FPlane::InitBounded()
{
  // L. Broglia

  projectedBoundary =
    surfaceBoundary.Project( pplace.GetToPlacementCoordinates() );
}

G4double G4FPlane::HowNear( const G4Vector3D& x ) const
{
  const G4Point3D Pt = x;
  //G4double d = ClosestDistanceToPoint(Pt);
  //return d;
  return ( Pt.x()*Pl.a + Pt.y()*Pl.b + Pt.z()*Pl.c - Pl.d);
}


