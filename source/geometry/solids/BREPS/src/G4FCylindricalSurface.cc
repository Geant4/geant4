// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FCylindricalSurface.cc,v 1.4 1999-01-20 07:34:49 broglia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/*  /usr/local/gismo/repo/geometry/FG4Cylinder.cc,v 1.1 1992/10/27 22:02:29 alanb Exp  */
//  File:  FG4Cylinder.cc
//  Author:  Alan Breakstone

//  Contents ----------------------------------------------------------
//
//      FG4Cylinder::FG4Cylinder( const G4Point3D& o, const G4ThreeVec& a,
//	                          G4double r, G4double l )
//	FG4Cylinder::FG4Cylinder( const FG4Cylinder& c )
//	FG4Cylinder::PrintOn( ostream& os ) const
//	FG4Cylinder::operator==( const FG4Cylinder& c )
//	FG4Cylinder::WithinBoundary( const G4ThreeVec& x ) const
//	FG4Cylinder::Scale() const
//	FG4Cylinder::resize( G4double r, G4double l )
//
//  End ---------------------------------------------------------------

#include "G4FCylindricalSurface.hh"
#include "G4Sort.hh"


G4FCylindricalSurface::G4FCylindricalSurface( const G4Point3D& o, 
					      const G4Vector3D& a,
					      const G4double r, 
					      const G4double l 
					    ) 
{ 
  //  make a G4FCylindricalSurface with origin o, axis a, 
  //  radius r, and length l
  G4Vector3D dir(1,1,1);
  Position.Init(dir, a, o);
  origin = o;

  
  //  Require length to be positive or zero
  //	if ( l > 0.0 )
  if ( l >= 0.0 )
    length = l;
  else 
  {
    G4cerr << "Error in G4FCylindricalSurface::G4FCylindricalSurface"
	   << "--asked for negative length\n"
	   << "\tDefault length of 0.0 is used.\n";

    length = 0.0;
  }

  //  Require radius to be non-negative (i.e., allow zero)
  if ( r >= 0.0 )
    radius = r;
  else 
  {
    G4cerr << "Error in G4FCylindricalSurface::G4FCylindricalSurface"
	   << "--asked for negative radius\n"
	   << "\tDefault value of 0.0 is used.\n";
    
    radius = 0.0;
  }

}


//  copy constructor
G4FCylindricalSurface::G4FCylindricalSurface( const G4FCylindricalSurface& c )
{ 
  length = c.length;
}

//  printing function using C++ ostream class
void G4FCylindricalSurface::PrintOn( ostream& os ) const
{ 
  os << "G4FCylindricalSurface with origin: " << origin << "\t"
     << "and axis: " << Position.GetAxis() << "\n"
     << "\t radius: " << radius << "\t and length: "
     << length << "\n";
}


int G4FCylindricalSurface::operator==( const G4FCylindricalSurface& c )
{
/*  return (   origin == c.origin && 
	     axis   == c.axis   && 
	     radius == c.radius && 
	     length == c.length    );*/
  return 1;
}


// Added 18.7-95
// Modified by L. Broglia (01/12/98)
void G4FCylindricalSurface::CalcBBox()
{
  // Finds the bounds of the surface iow
  // calculates the bounds for a bounding box
  // to the surface. The bounding box is used
  // for a preliminary check of intersection.
  G4Point3D Max = -PINFINITY;
  G4Point3D Min =  PINFINITY;

  G4Point3D Tmp; 
  G4Point3D Origin    = Position.GetLocation();
  G4Point3D EndOrigin = Origin + (length*Position.GetAxis());
  G4Point3D Radius(radius, radius, 0);
 
  // Default BBox
  G4Point3D Tolerance(kCarTolerance, kCarTolerance, kCarTolerance);
  G4Point3D BoxMin(Origin-Tolerance);
  G4Point3D BoxMax(Origin+Tolerance);

  bbox = new G4BoundingBox3D();
  bbox->Init(BoxMin, BoxMax);


  Tmp = (Origin - Radius);
  bbox->Extend(Tmp);
  
  Tmp = Origin + Radius;
  bbox->Extend(Tmp);

  Tmp = EndOrigin - Radius;
  bbox->Extend(Tmp);
  
  Tmp = EndOrigin + Radius;
  bbox->Extend(Tmp);
}


int G4FCylindricalSurface::Intersect( const G4Ray& ry )  
{
  //  Distance along a Ray (straight line with G4ThreeVec) to leave or enter
  //  a G4CylindricalSurface.  The input variable which_way should be set 
  //  to +1 to indicate leaving a G4CylindricalSurface, -1 to indicate 
  //  entering a G4CylindricalSurface.
  //  p is the point of intersection of the Ray with the G4CylindricalSurface.
  //  If the G4Vector3D of the Ray is opposite to that of the Normal to
  //  the G4CylindricalSurface at the intersection point, it will not leave 
  //  the G4CylindricalSurface.
  //  Similarly, if the G4Vector3D of the Ray is along that of the Normal 
  //  to the G4CylindricalSurface at the intersection point, it will not enter
  //  the G4CylindricalSurface.
  //  This method is called by all finite shapes sub-classed to 
  //  G4CylindricalSurface.
  //  Use the virtual function table to check if the intersection point
  //  is within the boundary of the finite shape.
  //  A negative result means no intersection.
  //  If no valid intersection point is found, set the distance
  //  and intersection point to large numbers.
  
  int which_way=1;
  
  if(!Inside(ry.GetStart()))
    which_way = -1;	   
  
  distance = FLT_MAXX;
  G4Vector3D lv ( FLT_MAXX, FLT_MAXX, FLT_MAXX );
  
  closest_hit = lv;

  //  Origin and G4Vector3D unit vector of Ray.
  G4Vector3D x    = ry.GetStart();
  G4Vector3D dhat = ry.GetDir();

  //  Axis unit vector of the G4CylindricalSurface.
  G4Vector3D ahat = GetAxis();
  
  //  array of solutions in distance along the Ray
  G4double s[2];
  s[0] = -1.0; 
  s[1] = -1.0 ;

  //  calculate the two solutions (quadratic equation)
  G4Vector3D d    = x - GetOrigin();
  G4double   radiu = GetRadius();

  //quit with no intersection if the radius of the G4CylindricalSurface is zero
  //	if ( radiu <= 0.0 )
  //		return 0;
	
  G4double dsq  = d * d;
  G4double da   = d * ahat;
  G4double dasq = da * da;
  G4double rsq  = radiu * radiu;
  G4double qsq  = dsq - dasq;
  G4double dira = dhat * ahat;
  G4double a    = 1.0 - dira * dira;
  
  if ( a <= 0.0 )
    return 0;
	
  G4double b       = 2. * ( d * dhat - da * dira );
  G4double c       = rsq - qsq;
  G4double radical = b * b + 4. * a * c; 
	
  if ( radical < 0.0 ) 
    return 0;

  G4double root = sqrt( radical );
  s[0] = ( - b + root ) / ( 2. * a );
  s[1] = ( - b - root ) / ( 2. * a );

  // Validity of the solutions
  // the hit point must be into the bounding box of the conical surface
  G4Point3D p0 = x + s[0]*dhat;
  G4Point3D p1 = x + s[1]*dhat;

  if( !GetBBox()->Inside(p0) )
    s[0] = kInfinity;

  if( !GetBBox()->Inside(p1) )
    s[1] = kInfinity;
  
  //  now loop over each positive solution, keeping the first one (smallest
  //  distance along the Ray) which is within the boundary of the sub-shape
  //  and which also has the correct G4Vector3D with respect to the Normal to
  //  the G4CylindricalSurface at the intersection point
  G4int nbinter = 0;
  distance = kInfinity;

  for ( G4int i = 0; i < 2; i++ ) 
  {  
    if(s[i] < kInfinity)
      if ( s[i] >= kCarTolerance*0.5 ) 
      {   
	if( distance > s[i]*s[i])
	  distance = s[i]*s[i];
	
	nbinter ++;
      }    
      else if ( s[i] >= -kCarTolerance*0.5 ) 
      {
	// the point is on the surface
	distance = 0;
	return 1;
      }   
  }

  return nbinter;
}


G4double G4FCylindricalSurface::HowNear( const G4Vector3D& x ) const
{
  //  Distance from the point x to the infinite G4CylindricalSurface.
  //  The distance will be positive if the point is Inside the 
  //  G4FCylindricalSurface, negative if the point is outside.

  G4Vector3D d = x - origin;
  G4double dA = d * Position.GetAxis();
  G4double rad = sqrt( d.mag2() - dA*dA );
  G4double hownear;

  // solve problem of tolerance
  if(fabs(dA) < kCarTolerance)
    dA = 0;
  
  if(dA > length)
    hownear = length - dA;
  else if(dA < 0)
    hownear = dA;
  else
    hownear = radius - rad;
 
  return hownear;
}

int G4FCylindricalSurface::WithinBoundary( const G4Vector3D& x ) const
{
  //  return 1 if point x is within the boundaries of the G4FCylindricalSurface
  //  return 0 otherwise (assume it is on the cylinder)
  if ( fabs( ( x - Position.GetLocation()) * Position.GetAxis() ) <= 0.5 * length )
    return 1;
  else
    return 0;
}


G4double G4FCylindricalSurface::Scale() const
{
  //  Returns the radius of a G4FCylindricalSurface unless it is zero, in which
  //  case returns the length.
  //  Used for Scale-invariant tests of surface thickness.
  if ( radius == 0.0 )
    return length;
  else
    return radius;
}


G4Vector3D G4FCylindricalSurface::SurfaceNormal( const G4Point3D& p ) const
{
  //  return the Normal unit vector to the G4CylindricalSurface at a point 
  //  p on (or nearly on) the G4CylindricalSurface
  
  G4Vector3D n = ( p - Position.GetLocation() ) - 
    ( ( p - Position.GetLocation()) * Position.GetAxis() ) *Position.GetAxis();
  G4double nmag = n.mag();
  
  if ( nmag != 0.0 )
    n = n * (1/nmag);
  
  return n;
}

int G4FCylindricalSurface::Inside ( const G4Vector3D& x ) const
{ 
  //  Return 0 if point x is outside G4CylindricalSurface, 1 if Inside.
  //  Outside means that the distance to the G4CylindricalSurface would 
  //  be negative.
  //  Use the HowNear function to calculate this distance.
  if ( HowNear( x ) >= -0.5*kCarTolerance )
    return 1;
  else
    return 0; 
}


void G4FCylindricalSurface::resize( G4double r, G4double l )
{
  //  Resize a G4FCylindricalSurface to a new radius r and new length l
  //  Require radius to be non-negative
  if ( r >= 0.0 )
    radius = r;
  else 
  {
    G4cerr << "Error in G4FCylindricalSurface::resize"
	   << "--asked for negative radius\n"
	   << "\tOriginal value of " << radius << " is retained.\n";
  }

  //  Require length to be positive
  if ( l > 0.0 )
    length = l;
  else 
  {
    G4cerr << "Error in G4FCylindricalSurface::resize"
	   << "--asked for negative or zero length\n"
	   << "\tOriginal value of " << length << " is retained.\n";
  }
}


