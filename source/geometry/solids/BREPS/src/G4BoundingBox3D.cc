// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BoundingBox3D.cc,v 1.5 2000-11-08 14:22:09 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4BoundingBox3D.cc
//
// ----------------------------------------------------------------------

#include "G4BoundingBox3D.hh"
#include "geomdefs.hh"

const G4BoundingBox3D G4BoundingBox3D::
          space( G4Point3D(-kInfinity, -kInfinity, -kInfinity),
		 G4Point3D(+kInfinity, +kInfinity, +kInfinity)  );

/////////////////////////////////////////////////////////////////////////////

G4BoundingBox3D::G4BoundingBox3D() { distance =0; }

G4BoundingBox3D::G4BoundingBox3D(const G4Point3D& p1, const G4Point3D& p2) 
{
  Init(p1, p2);
}

G4BoundingBox3D::G4BoundingBox3D(const G4Point3D& p) 
{
  Init(p);
}

G4BoundingBox3D::~G4BoundingBox3D() {}

G4BoundingBox3D::G4BoundingBox3D(const G4BoundingBox3D& right)
  : box_min(right.box_min), box_max(right.box_max),
    distance(right.distance), test_result(right.test_result),
    MiddlePoint(right.MiddlePoint), GeantBox(right.GeantBox)
{
}

G4BoundingBox3D& G4BoundingBox3D::operator=(const G4BoundingBox3D& right)
{
  if (&right == this) return *this;
  box_min  = right.box_min;
  box_max  = right.box_max;
  distance = right.distance;
  test_result = right.test_result;
  MiddlePoint = right.MiddlePoint;
  GeantBox = right.GeantBox;
  
  return *this;
}

void G4BoundingBox3D::Init(const G4Point3D& p1, const G4Point3D& p2) 
{
  // L. Broglia
  // Maybe temporary
  // Create a BBox bigger than the reality

  box_min.setX( G4std::min(p1.x(), p2.x()) - kCarTolerance );
  box_min.setY( G4std::min(p1.y(), p2.y()) - kCarTolerance );
  box_min.setZ( G4std::min(p1.z(), p2.z()) - kCarTolerance );
  box_max.setX( G4std::max(p1.x(), p2.x()) + kCarTolerance );
  box_max.setY( G4std::max(p1.y(), p2.y()) + kCarTolerance );
  box_max.setZ( G4std::max(p1.z(), p2.z()) + kCarTolerance );
  
  // Calc half spaces
  GeantBox = (box_max - box_min)*0.5;
  MiddlePoint = (box_min + box_max)*0.5;
  distance = 0;
}


void G4BoundingBox3D::Init(const G4Point3D& p) 
{
  box_min= box_max= MiddlePoint= p;
  GeantBox= G4Point3D(0, 0, 0);
  distance= 0;
}


/////////////////////////////////////////////////////////////////////////////

void G4BoundingBox3D::Extend(const G4Point3D& p) 
{
  
  // L. Broglia
  // Maybe temporary
  // Create a BBox bigger than the reality

  if (p.x() < box_min.x()) 
    box_min.setX( p.x() - kCarTolerance );
  else if (p.x() > box_max.x())
    box_max.setX( p.x() + kCarTolerance );
 
  if (p.y() < box_min.y()) 
    box_min.setY( p.y() - kCarTolerance );
  else if (p.y() > box_max.y()) 
    box_max.setY( p.y() + kCarTolerance );

  if (p.z() < box_min.z())
    box_min.setZ( p.z() - kCarTolerance );
  else if (p.z() > box_max.z())
    box_max.setZ( p.z() + kCarTolerance );

  // L. Broglia
  // Now re-calculate GeantBox and MiddlePoint
  GeantBox    = (box_max - box_min)*0.5;
  MiddlePoint = (box_min + box_max)*0.5;
  
}

////////////////////////////////////////////////////////////////////////////


G4int G4BoundingBox3D::Test(const G4Ray& rayref)
{
  const G4Point3D&  tmp_ray_start = rayref.GetStart();
  const G4Vector3D& tmp_ray_dir   = rayref.GetDir();

  G4Point3D  ray_start = tmp_ray_start ;
  G4Vector3D ray_dir   = tmp_ray_dir   ;

  G4double rayx,rayy,rayz;
  rayx = ray_start.x();
  rayy = ray_start.y();
  rayz = ray_start.z();

  // Test if ray starting point is in the bbox or not
  if((rayx < box_min.x()) || (rayx > box_max.x()) ||
     (rayy < box_min.y()) || (rayy > box_max.y()) ||		
     (rayz < box_min.z()) || (rayz > box_max.z())   )
  {
    // Outside, check for intersection with bbox
    
    // Adapt ray_starting point to box

    const G4Point3D ray_start2 = G4Point3D( ray_start - MiddlePoint );
    distance = DistanceToIn(ray_start2, ray_dir);

    if(!distance)
      test_result = 0; // Miss
    else
      test_result = 1; // Starting point outside box & hits box
  }
  else
  {
    // Inside
    // G4cout << "\nRay starting point Inside bbox.";
    test_result = 1;
    distance = 0;
  }

  return test_result;
}

///////////////////////////////////////////////////////////////////////////////


// Does an intersection exist?
//
// ALGORITHM:
//
// Check that if point lies outside x/y/z extent of box, travel is towards
// the box (ie. there is a possiblity of an intersection)


G4int G4BoundingBox3D::BoxIntersect(const G4Point3D&  gbox, 
				    const G4Point3D&  p   , 
				    const G4Vector3D& v    ) const
{
  G4double safx, safy, safz;
  G4double  fdx,  fdy,  fdz;
  
  fdx = GeantBox.x();    
  fdy = GeantBox.y();    
  fdz = GeantBox.z();

  safx=fabs(p.x())-fdx;   // minimum distance to x surface of shape
  safy=fabs(p.y())-fdy;
  safz=fabs(p.z())-fdz;
  
  // Will we Intersect?
  // If safx/y/z is >=0 the point is outside/on the box's x/y/z extent.
  // If both p.X()/y/z and v.X()/y/z repectively are both positive/negative,
  // travel is in a G4ThreeVec away from the shape.

  if ( ( (p.x()*v.x()>=0.0 ) && safx>0.0 ) || 
       ( (p.y()*v.y()>=0.0 ) && safy>0.0 ) ||
       ( (p.z()*v.z()>=0.0 ) && safz>0.0 )    )
    return 0; // No intersection  	
  else
    return 1; // Possible intersection
}

///////////////////////////////////////////////////////////////////////////////


// Distance to in
// Calculate distance to box from outside - return kBig if no intersection
//
// ALGORITHM:
//
// Check that if point lies outside x/y/z extent of box, travel is towards
// the box (ie. there is a possiblity of an intersection)
//
// Calculate pairs of minimum and maximum distances for x/y/z travel for
// intersection with the box's x/y/z extent.
// If there is a valid intersection, it is given by the maximum min distance
// (ie. distance to satisfy x/y/z intersections) *if* <= minimum max distance
// (ie. distance after which 1+ of x/y/z intersections not satisfied)
//
// NOTE:
//
// `Inside' safe - meaningful answers given if point is Inside the exact
// shape.

//G4double G4BoundingBox::distance_to_in(const G4Point3d& gbox, const G4Point3d& p, const G4ThreeVec& v) const
G4double G4BoundingBox3D::DistanceToIn(const G4Point3D& p, 
				       const G4Vector3D& v) const
{
    G4double safx,  safy,  safz,  snxt = 0;  // snxt = default return value
    G4double smin, sminx, sminy, sminz;
    G4double smax, smaxx, smaxy, smaxz;
    G4double stmp;
    G4double kBig = 10e20;
    G4double fdx,fdy,fdz;
    
    fdx = GeantBox.x();        
    fdy = GeantBox.y();        
    fdz = GeantBox.z();    

    safx = fabs(p.x())-fdx;   // minimum distance to x surface of shape
    safy = fabs(p.y())-fdy;
    safz = fabs(p.z())-fdz;

    // Will we Intersect?
    // If safx/y/z is >=0 the point is outside/on the box's x/y/z extent.
    // If both p.X()/y/z and v.X()/y/z repectively are both positive/negative,
    // travel is in a G4ThreeVec away from the shape.

    if ( ( ( p.x()*v.x()>=0.0 ) && safx>0.0) || 
	 ( ( p.y()*v.y()>=0.0 ) && safy>0.0) || 
	 ( ( p.z()*v.z()>=0.0 ) && safz>0.0)    )
      return snxt;   	
    
    // Compute min / max distance for x/y/z travel:
    if (safx<0.0)
    {
      // Inside x extent => Calc distance until trajectory leaves extent
      sminx=0.0;
      if (v.x()) 
	smaxx = fdx/fabs(v.x()) - p.x()/v.x();
      else
	smaxx = kBig;
    }
    else
    {
      // Outside extent or on boundary
      if (v.x()==0)
	return snxt; // Travel parallel
      else
      {
	stmp  = fabs(v.x());
	sminx = safx/stmp;
	smaxx = (fdx+fabs(p.x()))/stmp;
      }
    }
    
    if (safy<0.0)
    {
      // Inside y extent => Calc distance until trajectory leaves extent
      sminy=0.0;
      if (v.y()) 
	smaxy = fdy/fabs(v.y()) - p.y()/v.y();
      else
	smaxy = kBig;
    }
    else
    {
      // Outside extent or on boundary
      if (v.y()==0)
	return snxt; // Travel parallel
      else
      {
	stmp  = fabs(v.y());
	sminy = safy/stmp;
	smaxy = (fdy+fabs(p.y()))/stmp;
      }
    }
    
    if (safz<0.0)
    {
      // Inside z extent => Calc distance until trajectory leaves extent
      sminz=0.0;
      if (v.z()) 
	smaxz = fdz/fabs(v.z()) - p.z()/v.z();
      else 
	smaxz = kBig;
    }
    else
    {
      // Outside extent or on boundary
      if (v.z()==0)
	return snxt; // Travel parallel
      else
      {
	stmp  = fabs(v.z());
	sminz = safz/stmp;
	smaxz = (fdz+fabs(p.z()))/stmp;
      }
    }

    // Find minimum allowed Dist given min/max pairs
    if (sminx>sminy) 
      smin = sminx; // MAX(sminx,sminy,sminz)
    else 
      smin = sminy;
    
    if (sminz>smin) 
      smin=sminz;

    if (smaxx<smaxy) 
      smax = smaxx; // MIN(smaxx,smaxy,smaxz)
    else 
      smax = smaxy;
    
    if (smaxz<smax) 
      smax = smaxz;

    // If smin <= kCarTolerance then only clipping `tolerant' Area
    // -> no intersection
    
    G4double kCarTolerance = 0;
    
    if (smin>kCarTolerance && smin<=smax) 
      snxt=smin;
    
    return snxt;
}


///////////////////////////////////////////////////////////////////////////////

G4int G4BoundingBox3D::Inside(const G4Point3D& Pt) const
{
  if( ( Pt.x() >= box_min.x() && Pt.x() <= box_max.x() ) &&
      ( Pt.y() >= box_min.y() && Pt.y() <= box_max.y() ) &&
      ( Pt.z() >= box_min.z() && Pt.z() <= box_max.z() )    )
    return 1;
  else
    return 0;
}
