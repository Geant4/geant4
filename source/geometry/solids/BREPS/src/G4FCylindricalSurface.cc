// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FCylindricalSurface.cc,v 1.11 2000-11-08 14:22:10 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4FCylindricalSurface.cc
//
// ----------------------------------------------------------------------

#include "G4FCylindricalSurface.hh"
#include "G4Sort.hh"


G4FCylindricalSurface::G4FCylindricalSurface()
  : length(1.)
{
}


G4FCylindricalSurface::~G4FCylindricalSurface()
{
}


G4FCylindricalSurface::G4FCylindricalSurface( const G4Point3D& o, 
					      const G4Vector3D& a,
					      G4double r, 
					      G4double l 
					    ) 
{ 
  //  make a G4FCylindricalSurface with origin o, axis a, 
  //  radius r, and length l
  G4Vector3D dir(1,1,1);
  Position.Init(dir, a, o);

  origin = o;
  radius = r;
  
  //  Require length to be positive or zero
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


const char* G4FCylindricalSurface::NameOf() const 
{
  return "G4FCylindricalSurface"; 
}


void G4FCylindricalSurface::PrintOn( G4std::ostream& os ) const
{ 
  os << "G4FCylindricalSurface with origin: " << origin << "\t"
     << "and axis: " << Position.GetAxis() << "\n"
     << "\t radius: " << radius << "\t and length: "
     << length << "\n";
}


G4double G4FCylindricalSurface::Area() const 
{
  return ( 2.0 * M_PI * radius * length );
}


// Added 18.7-95
// Modified by L. Broglia (01/12/98)
void G4FCylindricalSurface::CalcBBox()
{
  // Finds the bounds of the surface iow
  // calculates the bounds for a bounding box
  // to the surface. The bounding box is used
  // for a preliminary check of intersection.
  G4Point3D Max = G4Point3D(-PINFINITY);
  G4Point3D Min = G4Point3D( PINFINITY);

  G4Point3D Tmp; 
  G4Point3D Origin    = Position.GetLocation();
  G4Point3D EndOrigin = G4Point3D( Origin + (length*Position.GetAxis()) );
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


G4int G4FCylindricalSurface::Intersect( const G4Ray& ry )  
{
  // This function count the number of intersections of a 
  // bounded cylindrical surface by a ray.
  // At first, calculates the intersections with the infinite 
  // cylindrical surfsace. After, count the intersections within the
  // finite cylindrical surface boundaries, and set "distance" to the 
  // closest distance from the start point to the nearest intersection
  // If the point is on the surface it returns or the intersection with
  // the opposite surface or kInfinity

  // If no intersection is founded, set distance = kInfinity and
  // return 0

  distance    = kInfinity;
  closest_hit = PINFINITY;

  // origin and direction of the ray
  G4Point3D  x    = ry.GetStart();
  G4Vector3D dhat = ry.GetDir();

  // cylinder axis 
  G4Vector3D ahat  = Position.GetAxis();
 
  //  array of solutions in distance along the ray
  G4double s[2];
  s[0]=-1.0;
  s[1]=-1.0;

  // calculate the two intersections (quadratic equation)   
  G4Vector3D gamma =  G4Vector3D( x - Position.GetLocation() );
  
  G4double ga = gamma * ahat;
  G4double da = dhat * ahat;

  G4double A = da * da - dhat * dhat;
  G4double B = 2 * ( -gamma * dhat + ga * da );
  G4double C =  -gamma * gamma + ga * ga  + radius * radius ;

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
  // the hit point must be into the bounding box of the cylindrical surface
  G4Point3D p0 = G4Point3D( x + s[0]*dhat );
  G4Point3D p1 = G4Point3D( x + s[1]*dhat );

  if( !GetBBox()->Inside(p0) )
    s[0] = kInfinity;

  if( !GetBBox()->Inside(p1) )
    s[1] = kInfinity;
  
  //  now loop over each positive solution, keeping the first one (smallest
  //  distance along the Ray) which is within the boundary of the sub-shape
  G4int nbinter = 0;
  distance = kInfinity;

  for ( G4int i = 0; i < 2; i++ ) 
  {  
    if(s[i] < kInfinity) {
      if ( s[i] >= kCarTolerance*0.5 ) {
	nbinter ++;
       	// real intersection
	// set the distance if it is the smallest
	if( distance > s[i]*s[i]) {
	  distance = s[i]*s[i];
	}
      }
    }    
  }

  return nbinter;
}


G4double G4FCylindricalSurface::HowNear( const G4Vector3D& x ) const
{
  // Shortest distance from the point x to the G4FCylindricalSurface.
  // The distance will be always positive 

  G4double   hownear;

  G4Vector3D upcorner = G4Vector3D ( radius, 0 , origin.z()+length);
  G4Vector3D downcorner = G4Vector3D ( radius, 0 , origin.z());
  G4Vector3D xd;  
  
  xd = G4Vector3D ( sqrt ( x.x()*x.x() + x.y()*x.y() ) , 0 , x.z() );
    
  
  G4double Zinter = (xd.z()) ;
  
  if ( ((Zinter >= downcorner.z()) && (Zinter <=upcorner.z())) ) {
    hownear = fabs( radius - xd.x() );
  } else {
    hownear = G4std::min ( (xd-upcorner).mag() , (xd-downcorner).mag() );
  }

  return hownear;
}

G4int G4FCylindricalSurface::WithinBoundary( const G4Vector3D& x ) const
{
  //  return 1 if point x is within the boundaries of the G4FCylindricalSurface
  //  return 0 otherwise (assume it is on the cylinder)
  if ( fabs( ( x - Position.GetLocation()) * Position.GetAxis() )
       <= 0.5 * length )
    return 1;
  else
    return 0;
}


G4double G4FCylindricalSurface::Scale() const
{
  //  Returns the radius of a G4FCylindricalSurface unless it is zero,
  //  in which case returns the length.
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
  
  G4Vector3D n = G4Vector3D( ( p - Position.GetLocation() ) - 
                           ( ( p - Position.GetLocation()) *
			       Position.GetAxis() ) * Position.GetAxis() );
  G4double nmag = n.mag();
  
  if ( nmag != 0.0 )
    n = n * (1/nmag);

  if( !sameSense )
    n = -n;
  
  return n;
}

G4int G4FCylindricalSurface::Inside ( const G4Vector3D& x ) const
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
