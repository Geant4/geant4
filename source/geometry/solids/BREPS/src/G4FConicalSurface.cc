// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FConicalSurface.cc,v 1.7 1999-01-27 16:12:02 broglia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/*  /usr/local/gismo/repo/geometry/G4FConicalSurface.cc,v 1.2 1993/02/05 00:38:39 alanb Exp  */
//  File:  G4FConicalSurface.cc
//  Author:  Alan Breakstone

//  Contents ----------------------------------------------------------
//
//      G4FConicalSurface::G4FConicalSurface( const G4Point3D& o, 
//                                            const G4Vector3D& a,
//		G4double l, G4double sr, G4double lr )
//	G4FConicalSurface::G4FConicalSurface( const G4FConicalSurface& c )
//	G4FConicalSurface::PrintOn( ostream& os ) const
//	G4FConicalSurface::operator==( const G4FConicalSurface& c )
//	G4FConicalSurface::WithinBoundary( const G4Vector3D& x ) const
//	G4FConicalSurface::Scale() const
//	G4FConicalSurface::Area() const
//	G4FConicalSurface::resize( G4double l, G4double sr, G4double lr )
//
//  End ---------------------------------------------------------------


#include "G4FConicalSurface.hh"
#include "G4Sort.hh"
#include "G4CircularCurve.hh"


G4FConicalSurface::G4FConicalSurface(const G4Point3D&  o, 
				     const G4Vector3D& a,
				     G4double          l, 
				     G4double          sr, 
				     G4double          lr 
				    ) 
{
  // Make a G4FConicalSurface with origin o, axis a, length l, small radius 
  // sr, and large radius lr. The angle is calculated below and the SetAngle
  // function of G4ConicalSurface is used to set it properly from the default
  // value used above in the initialization.
 
  // Create the position with origin o, axis a, and a direction
  G4Vector3D dir(1,1,1);
  Position.Init(dir, a, o);
  origin = o;
  
  //  Require length to be nonnegative
  if (l >=0)
    length = l;
  else 
  {
    G4cerr << "Error in G4FConicalSurface::G4FConicalSurface"
	   << "--asked for negative length\n"
	   << "\tDefault length of 0.0 is used.\n";
   
    length = 0.0;
  }
  
  //  Require small radius to be non-negative (i.e., allow zero)
  if ( sr >= 0.0 )
    small_radius = sr;
  else 
  {
    G4cerr << "Error in G4FConicalSurface::G4FConicalSurface"
	   << "--asked for negative small radius\n"
	   << "\tDefault value of 0.0 is used.\n";
   
    small_radius = 0.0;
  }

  //  Require large radius to exceed small radius
  if ( lr > small_radius )
    large_radius = lr;
  else 
  {
    G4cerr << "Error in G4FConicalSurface::G4FConicalSurface"
	   << "--large radius must exceed small radius\n"
	   << "\tDefault value of small radius +1 is used.\n";
   
    large_radius = small_radius + 1.0;
  }

  //  Calculate the angle of the G4ConicalSurface from the length and radii
  tan_angle =  ( large_radius - small_radius ) / length ;
}


G4FConicalSurface::G4FConicalSurface( const G4FConicalSurface& c )
  //: G4ConicalSurface( c.origin, c.GetAxis(), c.GetAngle() )
{ 
  //  copy constructor
  small_radius = c.small_radius;
  large_radius = c.large_radius;
  length       = c.length;
  tan_angle    = c.tan_angle;
}


// Modified by L. Broglia (01/12/98)
void G4FConicalSurface::CalcBBox()
{
  G4Point3D Max   = -PINFINITY;
  G4Point3D Min   =  PINFINITY;
  G4Point3D Tmp;
  G4double  delta = small_radius / tan_angle;

  G4Point3D Origin    = Position.GetLocation();
  G4Point3D EndOrigin = Origin + (length * Position.GetAxis());
  
  G4double radius = large_radius;
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


void G4FConicalSurface::PrintOn( ostream& os ) const
{ 
  //  printing function using C++ ostream class
  os << "G4FConicalSurface with origin: " << origin << "\t"
     << "and axis: " << Position.GetAxis() << "\n"
     << "\t small radius: " << small_radius 
     << "\t large radius: " << large_radius
     << "\t and length: " << length << "\n";
}


int G4FConicalSurface::operator==( const G4FConicalSurface& c )
{
  return ( origin             == c.origin                &&
	   Position.GetAxis() == c.Position.GetAxis()    &&
	   small_radius       == c.small_radius          && 
	   large_radius       == c.large_radius          && 
	   length      	      == c.length                &&
	   tan_angle          == c.tan_angle                );
}


int G4FConicalSurface::WithinBoundary( const G4Vector3D& x ) const
{ 
  //  return 1 if point x is within the boundaries of the G4FConicalSurface
  //  return 0 otherwise (assume it is on the G4ConicalSurface)
  G4Vector3D q = x - origin;
  
  G4double qmag = q.mag();
  G4double s    = sin( atan2(large_radius-small_radius, length) );
  G4double ls   = small_radius / s;
  G4double ll   = large_radius / s;
  
  if ( ( qmag >= ls )  &&  ( qmag <= ll ) )
    return 1;
  else
    return 0;
}


G4double G4FConicalSurface::Scale() const
{  
  //  Returns the small radius of a G4FConicalSurface unless it is zero, in 
  //  which case returns the large radius.
  //  Used for Scale-invariant tests of surface thickness.
  if ( small_radius == 0.0 )
    return large_radius;
  else
    return small_radius;
}


G4double G4FConicalSurface::Area() const
{ 
 //  Returns the Area of a G4FConicalSurface
  G4double rdif = large_radius - small_radius; 
  
  return ( M_PI * ( small_radius + large_radius ) * 
	   sqrt( length * length  +  rdif * rdif ) );
}


void G4FConicalSurface::resize( G4double l, G4double sr, G4double lr )
{
  //  Resize a G4FConicalSurface to a new length l, and new radii sr and lr.
  //  Must Reset angle of the G4ConicalSurface as well based on these new 
  //  values.
  //  Require length to be non-negative
  
  //	if ( l > 0.0 )
  if ( l >= 0.0 )
    length = l;
  else 
  {
    G4cerr << "Error in G4FConicalSurface::resize"
	   << "--asked for negative length\n"
	   << "\tOriginal value of " << length << " is retained.\n";
  }

  //  Require small radius to be non-negative (i.e., allow zero)
  if ( sr >= 0.0 )
    small_radius = sr;
  else 
  {
    G4cerr << "Error in G4FConicalSurface::resize"
	   << "--asked for negative small radius\n"
	   << "\tOriginal value of " << small_radius
	   << " is retained.\n";
  }

  //  Require large radius to exceed small radius
  if ( lr > small_radius )
    large_radius = lr;
  else 
  {
    G4double r = small_radius + 1.0;
    lr = ( large_radius <= small_radius ) ? r : large_radius;
    large_radius = lr;
    
    G4cerr << "Error in G4FConicalSurface::G4FConicalSurface"
	   << "--large radius must exceed small radius\n"
	   << "\tDefault value of " << large_radius << " is used.\n";
  }

  //  Calculate the angle of the G4ConicalSurface from the length and radii
  tan_angle =  ( large_radius - small_radius ) / length ;
 
}


int G4FConicalSurface::Intersect(const G4Ray& ry )
{ 
  // This function count the number of intersections of a 
  // bounded conical surface by a ray.
  // At first, calculates the intersections with the semi-infinite 
  // conical surfsace. After, count the intersections within the
  // finite conical surface boundaries, and set "distance" to the 
  // closest distance from the start point to the nearest intersection
  // If no intersection is founded, set distance = kInfinity and
  // return 0

  distance    = kInfinity;
  closest_hit = PINFINITY;

  // origin and direction of the ray
  G4Point3D  x    = ry.GetStart();
  G4Vector3D dhat = ry.GetDir();

  // cone angle and axis
  G4double   ta   = tan_angle;
  G4Vector3D ahat = Position.GetAxis();
 
  //  array of solutions in distance along the ray
  G4double s[2];
  s[0]=-1.0;
  s[1]=-1.0;

  // calculate the two intersections (quadratic equation)   
  G4Vector3D gamma =  x - Position.GetLocation();
  
  G4double t  = 1  +  ta * ta;
  G4double ga = gamma * ahat;
  G4double da = dhat * ahat;

  G4double A = t * da * da - dhat * dhat;
  G4double B = 2 * ( -gamma * dhat + t * ga * da - large_radius * ta * da);
  G4double C = ( -gamma * gamma + t * ga * ga 
		 - 2 * large_radius * ta * ga
		 + large_radius * large_radius );

  G4double radical = B * B  -  4.0 * A * C; 

  if ( radical < 0.0 ) 
    // no intersection
    return 0;
  else 
  {
    G4double root = sqrt( radical );
    s[0] = ( - B + root ) / ( 2. * A );
    s[1] = ( - B - root ) / ( 2. * A );
  }
  
  // validity of the solutions
  // the hit point must be into the bounding box of the conical surface
  G4Point3D p0 = x + s[0]*dhat;
  G4Point3D p1 = x + s[1]*dhat;
  
  if( !GetBBox()->Inside(p0) )
    s[0] = kInfinity;

  if( !GetBBox()->Inside(p1) )
    s[1] = kInfinity;
 
  // now loop over each positive solution, keeping the first one (smallest
  // distance along the ray) which is within the boundary of the sub-shape
  G4int nbinter = 0;
  distance = kInfinity;

  for ( G4int i = 0; i < 2; i++ ) 
  {  
    if(s[i] < kInfinity)
      if ( s[i] >= kCarTolerance*0.5 ) 
      { 
	// real intersection
	// set the distance if it is the smallest
	if( distance > s[i]*s[i])
	  distance = s[i]*s[i];
	
	// increase the number of real intersections
	nbinter ++;
      }    
      else if ( s[i] >= -kCarTolerance*0.5 ) 
      {
	// the point is on the surface
	// set the distance to 0 and return 1
	distance = 0;
	return 1;
      }   
  }

  return nbinter;
}


G4double G4FConicalSurface::HowNear( const G4Vector3D& x ) const
{ 
  // Shortest distance from the point x to the G4FConicalSurface.
  // The distance will be positive if the point is outside 
  // the G4ConicalSurface, negative if the point is inside.
  
  G4Vector3D d       = x - origin;
  G4double   dA      = d * Position.GetAxis();
  G4double   rad     = sqrt( d.mag2() - dA*dA );
  G4double   teta    = atan2( (large_radius - small_radius) , length );
  G4double   radiu   = rad - (large_radius - dA*tan_angle) ;
  G4double   hownear ;

  // solve problem of tolerance
  if(fabs(dA) < kCarTolerance)
    dA = 0;

  if (dA > length)
    // point outside, so hownear must be positive
    hownear = dA - length;
  else if (dA < 0)
    // point outside, so hownear must be positive
    hownear = -dA;
  else
    hownear = radiu * cos(teta);

  return hownear;
}


G4Vector3D G4FConicalSurface::SurfaceNormal( const G4Point3D& p ) const
{  
  //  return the Normal unit vector to the G4ConicalSurface at a point p 
  //  on (or nearly on) the G4ConicalSurface
  G4Vector3D s  = p - origin;
  G4double   da = s * Position.GetAxis();
  G4double   r  = sqrt( s*s - da*da);
  G4double   z  = tan_angle * r; 

  if (Position.GetAxis().z() < 0)
    z = -z; 

  G4Vector3D n(p.x(), p.y(), z);
  n = n.unit();
  
  if( !sameSense )
    n = -n;

  return n; 
}

int G4FConicalSurface::Inside ( const G4Vector3D& x ) const
{ 
  // Return 0 if point x is outside G4ConicalSurface, 1 if Inside.
  // Outside means that the distance to the G4ConicalSurface would be negative.
  // Use the HowNear function to calculate this distance.
  if ( HowNear( x ) >= -0.5*kCarTolerance )
    return 1;
  else
    return 0; 
}

