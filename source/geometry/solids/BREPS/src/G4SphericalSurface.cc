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
// $Id$
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4SphericalSurface.cc
//
// ----------------------------------------------------------------------

#include "G4SphericalSurface.hh"
#include "G4PhysicalConstants.hh"

G4SphericalSurface::G4SphericalSurface()
  : G4Surface(), radius(1.0), phi_1(0.0), phi_2(2*pi),
    theta_1(0.0), theta_2(pi)
{
   // default constructor
   // default x_axis is ( 1.0, 0.0, 0.0 ), z_axis is ( 0.0, 0.0, 1.0 ),
   // default radius is 1.0
   // default phi_1 is 0, phi_2 is 2*PI
   // default theta_1 is 0, theta_2 is PI

   x_axis = G4Vector3D( 1.0, 0.0, 0.0 );
   z_axis = G4Vector3D( 0.0, 0.0, 1.0 );
   //   OuterBoundary = new G4BREPPolyline();
}


G4SphericalSurface::G4SphericalSurface( const G4Vector3D&, 
                                        const G4Vector3D& xhat,
                                        const G4Vector3D& zhat,
                                        G4double r, 
                                        G4double ph1, G4double ph2, 
                                        G4double th1, G4double th2) 
  : G4Surface()
{ 
  // Require both x_axis and z_axis to be unit vectors

  G4double xhatmag = xhat.mag();
  if ( xhatmag != 0.0 )
  {
    x_axis = xhat * (1/ xhatmag); // this makes the x_axis a unit vector
  }
  else
  {
    std::ostringstream message;
    message << "x_axis has zero length." << G4endl
	    << "Default x_axis of (1, 0, 0) is used.";    
    G4Exception("G4SphericalSurface::G4SphericalSurface()",
                "GeomSolids1001", JustWarning, message);

    x_axis = G4Vector3D( 1.0, 0.0, 0.0 );
  }

  G4double zhatmag = zhat.mag();
  
  if (zhatmag != 0.0)
  {
    z_axis = zhat *(1/ zhatmag);  // this makes the z_axis a unit vector
  }
  else 
  {
    std::ostringstream message;
    message << "z_axis has zero length." << G4endl
	    << "Default z_axis of (0, 0, 1) is used.";    
    G4Exception("G4SphericalSurface::G4SphericalSurface()",
                "GeomSolids1001", JustWarning, message);

    z_axis = G4Vector3D( 0.0, 0.0, 1.0 );
  }
  
  //  Require radius to be non-negative
  //
  if ( r >= 0.0 )
  {
    radius = r;
  }
  else 
  {
    std::ostringstream message;
    message << "Radius cannot be less than zero." << G4endl
	    << "Default radius of 1.0 is used.";    
    G4Exception("G4SphericalSurface::G4SphericalSurface()",
                "GeomSolids1001", JustWarning, message);

    radius = 1.0;
  }

  //  Require phi_1 in the range: 0 <= phi_1 < 2*PI
  //  and phi_2 in the range: phi_1 < phi_2 <= phi_1 + 2*PI
  //
  if ( ( ph1 >= 0.0 ) && ( ph1 < 2*pi ) )
  {
    phi_1 = ph1;
  }
  else 
  {
    std::ostringstream message;
    message << "Lower azimuthal limit is out of range." << G4endl
	    << "Default angle of 0 is used.";    
    G4Exception("G4SphericalSurface::G4SphericalSurface()",
                "GeomSolids1001", JustWarning, message);

    phi_1 = 0.0;
  }

  if ( ( ph2 > phi_1 ) && ( ph2 <=  ( phi_1 + twopi ) ) )
  {
    phi_2 = ph2;
  }
  else 
  {
    std::ostringstream message;
    message << "Upper azimuthal limit is out of range." << G4endl
	    << "Default angle of 2*PI is used.";    
    G4Exception("G4SphericalSurface::G4SphericalSurface()",
                "GeomSolids1001", JustWarning, message);

    phi_2 = twopi;
  }
 
  //  Require theta_1 in the range: 0 <= theta_1 < PI 
  //  and theta-2 in the range: theta_1 < theta_2 <= theta_1 + PI
  //
  if ( ( th1 >= 0.0 ) && ( th1 < pi ) )
  {
    theta_1 = th1;
  }
  else
  {
    std::ostringstream message;
    message << "Lower polar limit is out of range." << G4endl
	    << "Default angle of 0 is used.";    
    G4Exception("G4SphericalSurface::G4SphericalSurface()",
                "GeomSolids1001", JustWarning, message);

    theta_1 = 0.0;
  }  
  
  if  ( ( th2 > theta_1 ) && ( th2 <= ( theta_1 + pi ) ) )
  {
    theta_2 =th2;
  }
  else 
  {
    std::ostringstream message;
    message << "Upper polar limit is out of range." << G4endl
	    << "Default angle of PI is used.";    
    G4Exception("G4SphericalSurface::G4SphericalSurface()",
                "GeomSolids1001", JustWarning, message);

    theta_2 = pi;
  } 
}


G4SphericalSurface::~G4SphericalSurface()
{
}


G4SphericalSurface::G4SphericalSurface( const G4SphericalSurface& surf )
  : G4Surface()
{
  x_axis = surf.x_axis;
  z_axis = surf.z_axis;
  radius = surf.radius;
  phi_1 = surf.phi_1;
  phi_2 = surf.phi_2;
  theta_1 = surf.theta_1;
  theta_2 = surf.theta_2;
}                               

G4SphericalSurface&
G4SphericalSurface::operator=( const G4SphericalSurface& surf )
{
  if (&surf == this)  { return *this; }
  x_axis = surf.x_axis;
  z_axis = surf.z_axis;
  radius = surf.radius;
  phi_1 = surf.phi_1;
  phi_2 = surf.phi_2;
  theta_1 = surf.theta_1;
  theta_2 = surf.theta_2;

  return *this;
}

const char* G4SphericalSurface::NameOf() const
{
  return "G4SphericalSurface";
}

void G4SphericalSurface::PrintOn( std::ostream& os ) const
{  
  os << "G4SphericalSurface surface with origin: " << origin << "\t"
     << "radius: " << radius << "\n"
     << "\t local x_axis: " << x_axis
     << "\t local z_axis: " << z_axis << "\n" 
     << "\t lower azimuthal limit: " << phi_1 << " radians\n"
     << "\t upper azimuthal limit: " << phi_2 << " radians\n"
     << "\t lower polar limit    : " << theta_1 << " radians\n"
     << "\t upper polar limit    : " << theta_2 << " radians\n";
}


G4double G4SphericalSurface::HowNear( const G4Vector3D& x ) const
{
  //  Distance from the point x to the G4SphericalSurface.
  //  The distance will be positive if the point is Inside the 
  //  G4SphericalSurface, negative if the point is outside.

  G4Vector3D d = G4Vector3D( x - origin );
  G4double rds = d.mag();

  return (radius - rds);
}


/*
G4double G4SphericalSurface::distanceAlongRay( G4int which_way,
                                               const G4Ray* ry,
                                               G4Vector3D& p ) const

 // Distance along a Ray (straight line with G4Vector3D) to leave or enter
 // a G4SphericalSurface.  The input variable which_way should be set to +1 to
 // indicate leaving a G4SphericalSurface, -1 to indicate entering the surface.
 // p is the point of intersection of the Ray with the G4SphericalSurface.
 // If the G4Vector3D of the Ray is opposite to that of the Normal to
 // the G4SphericalSurface at the intersection point, it will not leave the
 // G4SphericalSurface.
 // Similarly, if the G4Vector3D of the Ray is along that of the Normal 
 // to the G4SphericalSurface at the intersection point, it will not enter the
 // G4SphericalSurface.
 // This method is called by all finite shapes sub-classed to
 // G4SphericalSurface.
 // Use the virtual function table to check if the intersection point
 // is within the boundary of the finite shape.
 // A negative result means no intersection.
 // If no valid intersection point is found, set the distance
 // and intersection point to large numbers.

{
  G4double Dist = kInfinity;
  G4Vector3D lv ( kInfinity, kInfinity, kInfinity );
  p = lv;

  // Origin and G4Vector3D unit vector of Ray.
  //
  G4Vector3D x = ry->Position( 0.0 );
  G4Vector3D dhat = ry->Direction( 0.0 );
  G4int isoln = 0, maxsoln = 2;

  // Array of solutions in distance along the Ray
  //
  G4double s[2];s[0] = -1.0; s[1]= -1.0 ;

  // Calculate the two solutions (quadratic equation)
  //
  G4Vector3D d = x - GetOrigin();
  G4double radius = GetRadius();

  // Quit with no intersection if the radius of the G4SphericalSurface is zero
  //
  if ( radius <= 0.0 )  { return Dist; }

  G4double dsq = d * d;
  G4double rsq = radius * radius;
  G4double b = d * dhat;
  G4double c = dsq - rsq;
  G4double radical = b * b - c;

  // Quit with no intersection if the radical is negative
  //
  if ( radical < 0.0 )  { return Dist; }

  G4double root = std::sqrt( radical );
  s[0] = -b + root;
  s[1] = -b - root;

  // Order the possible solutions by increasing distance along the Ray
  //
  G4Sort_double( s, isoln, maxsoln-1 );

  // Now loop over each positive solution, keeping the first one (smallest
  // distance along the Ray) which is within the boundary of the sub-shape
  // and which also has the correct G4Vector3D with respect to the Normal to
  // the G4SphericalSurface at the intersection point
  //
  for ( isoln = 0; isoln < maxsoln; isoln++ )
  {
    if ( s[isoln] >= 0.0 )
    {
      if ( s[isoln] >= kInfinity )  { return Dist; }  // quit if too large
      Dist = s[isoln];
      p = ry->Position( Dist );
      if ( ( ( dhat * Normal( p ) * which_way ) >= 0.0 )
          && ( WithinBoundary( p ) == 1 ) )  { return Dist; }
    }
  }

  // Get here only if there was no solution within the boundary,
  // reset distance and intersection point to large numbers

  p = lv;
  return kInfinity;
}          
*/


void G4SphericalSurface::CalcBBox()
{
  G4double x_min = origin.x() - radius;
  G4double y_min = origin.y() - radius;
  G4double z_min = origin.z() - radius;
  G4double x_max = origin.x() + radius;
  G4double y_max = origin.y() + radius;
  G4double z_max = origin.z() + radius;  
    
  G4Point3D Min(x_min, y_min, z_min);
  G4Point3D Max(x_max, y_max, z_max);  
  bbox = new G4BoundingBox3D( Min, Max); 
}


G4int G4SphericalSurface::Intersect( const G4Ray& ry )
{
  //  Distance along a Ray (straight line with G4Vector3D) to leave or enter
  //  a G4SphericalSurface.  The input variable which_way should be set to +1 
  //  to indicate leaving a G4SphericalSurface, -1 to indicate entering a 
  //  G4SphericalSurface.
  //  p is the point of intersection of the Ray with the G4SphericalSurface.
  //  If the G4Vector3D of the Ray is opposite to that of the Normal to
  //  the G4SphericalSurface at the intersection point, it will not leave the
  //  G4SphericalSurface.
  //  Similarly, if the G4Vector3D of the Ray is along that of the Normal 
  //  to the G4SphericalSurface at the intersection point, it will not enter 
  //  the G4SphericalSurface.
  //  This method is called by all finite shapes sub-classed to 
  //  G4SphericalSurface.
  //  Use the virtual function table to check if the intersection point
  //  is within the boundary of the finite shape.
  //  A negative result means no intersection.
  //  If no valid intersection point is found, set the distance
  //  and intersection point to large numbers.

  G4int which_way = (G4int)HowNear(ry.GetStart());
  
  if(!which_way)  { which_way =-1; }
  distance = kInfinity;
  
  G4Vector3D lv ( kInfinity, kInfinity, kInfinity );

  // p = lv;
  //
  closest_hit = lv;

  // Origin and G4Vector3D unit vector of Ray.
  //
  G4Vector3D x= G4Vector3D( ry.GetStart() );
  G4Vector3D dhat = ry.GetDir();
  G4int isoln = 0, maxsoln = 2;

  // Array of solutions in distance along the Ray
  //
  G4double ss[2];
  ss[0] = -1.0 ; 
  ss[1] = -1.0 ;

  // Calculate the two solutions (quadratic equation)
  //
  G4Vector3D d = G4Vector3D( x - GetOrigin() );
  G4double r = GetRadius();

  // Quit with no intersection if the radius of the G4SphericalSurface is zero
  //
  if ( r <= 0.0 )  { return 0; }
  
  G4double dsq     = d * d;
  G4double rsq     = r * r;
  G4double b       = d * dhat;
  G4double c       = dsq - rsq;
  G4double radical = b * b - c;
  
  // Quit with no intersection if the radical is negative
  //
  if ( radical < 0.0 )  { return 0; }
  
  G4double root = std::sqrt( radical );
  ss[0] = -b + root;
  ss[1] = -b - root;

  // Order the possible solutions by increasing distance along the Ray
  // G4Sort_double( ss, isoln, maxsoln-1 );
  //        
  if(ss[0] > ss[1])
  {
    G4double tmp =ss[0];
    ss[0] = ss[1];
    ss[1] = tmp;
  }

  // Now loop over each positive solution, keeping the first one (smallest
  // distance along the Ray) which is within the boundary of the sub-shape
  // and which also has the correct G4Vector3D with respect to the Normal to
  // the G4SphericalSurface at the intersection point
  //
  for ( isoln = 0; isoln < maxsoln; isoln++ ) 
  {
    if ( ss[isoln] >= kCarTolerance*0.5 ) 
    {
      if ( ss[isoln] >= kInfinity )  { return 0; }  // quit if too large
      distance = ss[isoln];
      closest_hit = ry.GetPoint( distance );
      if ( ( ( dhat * Normal( closest_hit ) * which_way ) >= 0.0 ) && 
           ( WithinBoundary( closest_hit ) == 1 )                      )
      {
        distance =  distance*distance;    
        return 1;
      }
    }
  }
  
  // Get here only if there was no solution within the boundary,
  // reset distance and intersection point to large numbers
  //
  distance = kInfinity;
  closest_hit = lv;

  return 0;
}          


/*
G4double G4SphericalSurface::distanceAlongHelix( G4int which_way,
                                                 const Helix* hx,
                                                 G4Vector3D& p ) const
{
  // Distance along a Helix to leave or enter a G4SphericalSurface.
  // The input variable which_way should be set to +1 to
  // indicate leaving a G4SphericalSurface, -1 to indicate entering a
  // G4SphericalSurface.
  // p is the point of intersection of the Helix with the G4SphericalSurface.
  // If the G4Vector3D of the Helix is opposite to that of the Normal to
  // the G4SphericalSurface at the intersection point, it will not leave the
  // G4SphericalSurface.
  // Similarly, if the G4Vector3D of the Helix is along that of the Normal 
  // to the G4SphericalSurface at the intersection point, it will not enter
  // the G4SphericalSurface.
  // This method is called by all finite shapes sub-classed to
  // G4SphericalSurface.
  // Use the virtual function table to check if the intersection point
  // is within the boundary of the finite shape.
  // If no valid intersection point is found, set the distance
  // and intersection point to large numbers.
  // Possible negative distance solutions are discarded.

{
  G4double Dist = kInfinity;
  G4Vector3D lv ( kInfinity, kInfinity, kInfinity );
  p = lv;
  G4int isoln = 0, maxsoln = 4;

  // Array of solutions in turning angle
  //
  G4double s[4];s[0] = -1.0; s[1]= -1.0 ;s[2] = -1.0; s[3]= -1.0 ;

  // Helix parameters
  //
  G4double rh = hx->GetRadius();         // radius of Helix
  G4Vector3D oh = hx->position( 0.0 );   // origin of Helix
  G4Vector3D dh = hx->direction( 0.0 );  // initial G4Vector3D of Helix
  G4Vector3D prp = hx->getPerp();        // perpendicular vector
  G4double prpmag = prp.mag();
  G4double rhp = rh / prpmag;
  G4double rs = GetRadius();             // radius of G4SphericalSurface
  if ( rs == 0.0 )   { return Dist; }    // quit if zero radius
  G4Vector3D os = GetOrigin();           // origin of G4SphericalSurface

  // Calculate quantities of use later on
  //
  G4Vector3D alpha = rhp * prp;
  G4Vector3D beta = rhp * dh;
  G4Vector3D gamma = oh - os;

  // Only consider approximate solutions to quadratic order in the turning
  // angle along the Helix
  //
  G4double A = beta * beta  +  gamma * alpha;
  G4double B = 2.0 * gamma * beta;
  G4double C = gamma * gamma  -  rs * rs;

  // Case if quadratic term is zero
  //
  if ( std::fabs( A ) < kCarTolerance )
  {
    if ( B == 0.0 )  { return Dist; }      // no intersection, quit
    else             { s[0] = -C / B; }
  }
  else   // General quadratic solution, A != 0
  {
    G4double radical = B * B  -  4.0 * A * C;
    if ( radical < 0.0 )  { return Dist; } // no intersection, quit
                       
    G4double root = std::sqrt( radical );
    s[0] = ( -B + root ) / ( 2.0 * A );
    s[1] = ( -B - root ) / ( 2.0 * A );     
    if ( rh < 0.0 )
    {
      s[0] = -s[0];
      s[1] = -s[1];
    }
    s[2] = s[0] + twopi;
    s[3] = s[1] + twopi;
  }

  // Order the possible solutions by increasing turning angle
  //
  G4Sort_double( s, isoln, maxsoln-1 );

  // Now loop over each positive solution, keeping the first one (smallest
  // distance along the Helix) which is within the boundary of the sub-shape.
  //
  for ( isoln = 0; isoln < maxsoln; isoln++ )
  {
    if ( s[isoln] >= 0.0 )
    {
      // Calculate distance along Helix and position and G4Vector3D vectors.
      //
      Dist = s[isoln] * std::fabs( rhp );
      p = hx->position( Dist );
      G4Vector3D d = hx->direction( Dist );

      // Now do approximation to get remaining distance to correct this
      // solution iterate it until the accuracy is below the user-set
      // surface precision
      //
      G4double delta = 0.;  
      G4double delta0 = kInfinity;
      G4int dummy = 1;
      G4int iter = 0;
      G4int in0 = Inside( hx->position ( 0.0 ) );
      G4int in1 = Inside( p );
      G4double sc = Scale();
      while ( dummy )
      {
        iter++;

        // Terminate loop after 50 iterations and Reset distance to large
        // number, indicating no intersection with G4SphericalSurface. This
        // generally occurs if the Helix curls too tightly to intersect it.
        //
        if ( iter > 50 )
        {
          Dist = kInfinity;
          p = lv;
          break;
        }

        // Find distance from the current point along the above-calculated
        // G4Vector3D using a Ray.
        // The G4Vector3D of the Ray and the Sign of the distance are
        // determined by whether the starting point of the Helix is Inside
        // or outside of the G4SphericalSurface.
        //
        in1 = Inside( p );
        if ( in1 )    //  current point Inside
        {
          if ( in0 )    //  starting point Inside
          {
            Ray* r = new Ray( p, d );
            delta = distanceAlongRay( 1, r, p );
            delete r;
          }
          else         //  starting point outside
          {
            Ray* r = new Ray( p, -d );
            delta = -distanceAlongRay( 1, r, p );
            delete r;
          }
        }
        else          //  current point outside
        {
          if ( in0 )    //  starting point Inside
          {
            Ray* r = new Ray( p, -d );
            delta = -distanceAlongRay( -1, r, p );
            delete r;
          }
          else          //  starting point outside
          {
            Ray* r = new Ray( p, d );
            delta = distanceAlongRay( -1, r, p );
            delete r;
          }
        }

        // Test if distance is less than the surface precision, if so
        // terminate loop.
        //
        if ( std::fabs( delta / sc ) <= SURFACE_PRECISION )  { break; }

        // If delta has not changed sufficiently from the previous iteration, 
        // skip out of this loop.
        //
        if (std::fabs( ( delta-delta0 )/sc ) <= SURFACE_PRECISION)  { break; }

        // If delta has increased in absolute value from the previous iteration
        // either the Helix doesn't Intersect the G4SphericalSurface or the
        // approximate solution is too far from the real solution.  Try groping
        // for a solution. If not found, Reset distance to large number,
        // indicating no intersection with the G4SphericalSurface.
        //
        if ( ( std::fabs( delta ) > std::fabs( delta0 ) ) )
        {
          Dist = std::fabs( rhp ) * gropeAlongHelix( hx );
          if ( Dist < 0.0 )
          {
            Dist = kInfinity;
            p = lv;
          }
          else
          {
            p = hx->position( Dist );
          }
          break;
        }

        // Set old delta to new one.
        //
        delta0 = delta;

        // Add distance to G4SphericalSurface to distance along Helix.
        //
        Dist += delta;

        // Negative distance along Helix means Helix doesn't Intersect
        // G4SphericalSurface. Reset distance to large number, indicating
        // no intersection with G4SphericalSurface.
        //
        if ( Dist < 0.0 )
        {
          Dist = kInfinity;
          p = lv;
          break;
        }

        // Recalculate point along Helix and the G4Vector3D.
        //
        p = hx->position( Dist );
        d = hx->direction( Dist );
      }  //  end of while loop

      // Now have best value of distance along Helix and position for this
      // solution, so test if it is within the boundary of the sub-shape
      // and require that it point in the correct G4Vector3D with respect to
      // the Normal to the G4SphericalSurface.

      if ( ( Dist < kInfinity ) &&
           ( ( hx->direction( Dist ) * Normal( p ) * which_way ) >= 0.0 )
          && ( WithinBoundary( p ) == 1 ) )  { return Dist; }
    }  // end of if s[isoln] >= 0.0 condition
  }  // end of for loop over solutions

  // If one gets here, there is no solution, so set distance along Helix
  // and position to large numbers.
  //
  Dist = kInfinity;
  p = lv;

  return Dist;
}


G4Vector3D G4SphericalSurface::Normal( const G4Vector3D& p ) const
{
  //  Return the Normal unit vector to the G4SphericalSurface at a point p on
  //  (or nearly on) the G4SphericalSurface.

  G4Vector3D n = p - origin;
  G4double nmag = n.mag();

  //  If the point p happens to coincide with the origin (which is possible
  //  if the radius is zero), set the Normal to the z-axis unit vector.
  //
  if ( nmag != 0.0 )  { n = n / nmag; }
  else  { n = G4Vector3D( 0.0, 0.0, 1.0 ); }

  return n;
}
*/

         
G4Vector3D G4SphericalSurface::Normal( const G4Vector3D& p ) const
{ 
  // Return the Normal unit vector to the G4SphericalSurface at a point p on
  // (or nearly on) the G4SphericalSurface.

  G4Vector3D n = G4Vector3D( p - origin );
  G4double nmag = n.mag();
  
  //  If the point p happens to coincide with the origin (which is possible
  //  if the radius is zero), set the Normal to the z-axis unit vector.
  //
  if ( nmag != 0.0 )  { n = n * (1/ nmag); }
  else  { n = G4Vector3D( 0.0, 0.0, 1.0 ); }

  return n;
}


G4Vector3D G4SphericalSurface::SurfaceNormal( const G4Point3D& p ) const
{ 
  // Return the Normal unit vector to the G4SphericalSurface at a point p on
  // (or nearly on) the G4SphericalSurface.

  G4Vector3D n = G4Vector3D( p - origin );
  G4double nmag = n.mag();
  
  //  If the point p happens to coincide with the origin (which is possible
  //  if the radius is zero), set the Normal to the z-axis unit vector.
  //
  if ( nmag != 0.0 )  { n = n * (1/ nmag); }
  else  { n = G4Vector3D( 0.0, 0.0, 1.0 ); }

  return n;
}


G4int G4SphericalSurface::Inside ( const G4Vector3D& x ) const
{
  // Return 0 if point x is outside G4SphericalSurface, 1 if Inside.
  // Outside means that the distance to the G4SphericalSurface would 
  // be negative.
  // Use the HowNear function to calculate this distance.

  if ( HowNear( x ) >= 0.0 )  { return 1; }
  else                        { return 0; }
}


G4int G4SphericalSurface::WithinBoundary( const G4Vector3D& x ) const
{ 
  // Return 1 if point x is on the G4SphericalSurface, otherwise return zero
  // (x is assumed to lie on the surface of the G4SphericalSurface, so one 
  // only checks the angular limits)

  G4Vector3D y_axis = G4Vector3D( z_axis.cross( x_axis ) );

  // Components of x in the local coordinate system of the G4SphericalSurface
  //
  G4double px = x * x_axis;
  G4double py = x * y_axis;
  G4double pz = x * z_axis;

  // Check if within polar angle limits
  //
  G4double theta = std::acos( pz / x.mag() );  // acos in range 0 to PI
  
  // Normal case
  //
  if ( theta_2 <= pi ) 
  {
    if ( ( theta < theta_1 ) || ( theta > theta_2 ) )  { return 0; }
  }
  else 
  {
    theta += pi;
    if ( ( theta < theta_1 ) || ( theta > theta_2 ) )  { return 0; }
  }

  // Now check if within azimuthal angle limits
  //
  G4double phi = std::atan2( py, px );  // atan2 in range -PI to PI
  
  if ( phi < 0.0 )  { phi += twopi; }
  
  // Normal case
  //
  if ( ( phi >= phi_1 )  &&  ( phi <= phi_2 ) )  { return 1; }
  
  // This is for the case that phi_2 is greater than 2*PI
  //
  phi += twopi;
  
  if ( ( phi >= phi_1 )  &&  ( phi <= phi_2 ) )  { return 1; }

  // Get here if not within azimuthal limits
 
  return 0;
}


G4double G4SphericalSurface::Scale() const
{
  // Returns the radius of a G4SphericalSurface unless it is zero, in which
  // case returns the arbitrary number 1.0.
  // Used for Scale-invariant tests of surface thickness.

  if ( radius == 0.0 )  { return 1.0; }
  else                  { return radius; }
}


G4double G4SphericalSurface::Area() const
{
  // Returns the Area of a G4SphericalSurface

  return ( 2.0*( theta_2 - theta_1 )*( phi_2 - phi_1)*radius*radius/pi );
}


void G4SphericalSurface::resize( G4double r, 
                                 G4double ph1, G4double ph2,
                                 G4double th1, G4double th2 )
{ 
  // Resize the G4SphericalSurface to new radius r, new lower and upper 
  // azimuthal angle limits ph1 and ph2, and new lower and upper polar 
  // angle limits th1 and th2.
  
  //  Require radius to be non-negative
  //
  if ( r >= 0.0 )  { radius = r; }
  else 
  {
    std::ostringstream message;
    message << "Radius cannot be less than zero." << G4endl
	    << "Original value of " << radius << " is retained.";    
    G4Exception("G4SphericalSurface::resize()",
                "GeomSolids1001", JustWarning, message);
  }

  //  Require azimuthal angles to be within bounds
  
  if ( ( ph1 >= 0.0 ) && ( ph1 < twopi ) )  { phi_1 = ph1; }
  else 
  {
    std::ostringstream message;
    message << "Lower azimuthal limit out of range." << G4endl
	    << "Original value of " << phi_1 << " is retained.";    
    G4Exception("G4SphericalSurface::resize()",
                "GeomSolids1001", JustWarning, message);
  }
  
  if ( ( ph2 > phi_1 ) && ( ph2 <= ( phi_1 + twopi ) ) )  { phi_2 = ph2; }
  else 
  {
    ph2 = ( phi_2 <= phi_1 ) ? ( phi_1 + kAngTolerance ) : phi_2;
    phi_2 = ph2;
    std::ostringstream message;
    message << "Upper azimuthal limit out of range." << G4endl
	    << "Value of " << phi_2 << " is used.";    
    G4Exception("G4SphericalSurface::resize()",
                "GeomSolids1001", JustWarning, message);
  }

  // Require polar angles to be within bounds
  //
  if ( ( th1 >= 0.0 ) && ( th1 < pi ) )  { theta_1 = th1; }
  else 
  {
    std::ostringstream message;
    message << "Lower polar limit out of range." << G4endl
	    << "Original value of " << theta_1 << " is retained.";    
    G4Exception("G4SphericalSurface::resize()",
                "GeomSolids1001", JustWarning, message);
  }
  
  if ( ( th2 > theta_1 ) && ( th2 <= ( theta_1 + pi ) ) )  { theta_2 = th2; }
  else
  {
    th2 = ( theta_2 <= theta_1 ) ? ( theta_1 + kAngTolerance ) : theta_2;
    theta_2 = th2;
    std::ostringstream message;
    message << "Upper polar limit out of range." << G4endl
	    << "Value of " << theta_2 << " is used.";    
    G4Exception("G4SphericalSurface::resize()",
                "GeomSolids1001", JustWarning, message);
  }
}


/*
void G4SphericalSurface::
rotate( G4double alpha, G4double beta,
        G4double gamma, G4ThreeMat& m, G4int inverse )
{
  // Rotate G4SphericalSurface first about global x_axis by angle alpha,
  //                  second about global y-axis by angle beta,
  //               and third about global z_axis by angle gamma
  // by creating and using G4ThreeMat objects in Surface::rotate
  // angles are assumed to be given in radians
  // if inverse is non-zero, the order of rotations is reversed
  // the axis is rotated here, the origin is rotated by calling
  // Surface::rotate

  G4Surface::rotate( alpha, beta, gamma, m, inverse );
  x_axis = m * x_axis;
  z_axis = m * z_axis;
}


void G4SphericalSurface::
rotate( G4double alpha, G4double beta, G4double gamma, G4int inverse )
{
  // Rotate G4SphericalSurface first about global x_axis by angle alpha,
  //                  second about global y-axis by angle beta,
  //               and third about global z_axis by angle gamma
  // by creating and using G4ThreeMat objects in Surface::rotate
  // angles are assumed to be given in radians
  // if inverse is non-zero, the order of rotations is reversed
  // the axis is rotated here, the origin is rotated by calling
  // Surface::rotate

  G4ThreeMat m;
  G4Surface::rotate( alpha, beta, gamma, m, inverse );
  x_axis = m * x_axis;
  z_axis = m * z_axis;
}
*/


/*
G4double G4SphericalSurface::gropeAlongHelix( const Helix* hx ) const
{
  // Grope for a solution of a Helix intersecting a G4SphericalSurface.
  // This function returns the turning angle (in radians) where the
  // intersection occurs with only positive values allowed, or -1.0 if
  // no intersection is found.
  // The idea is to start at the beginning of the Helix, then take steps
  // of some fraction of a turn.  If at the end of a Step, the current position
  // along the Helix and the previous position are on opposite sides of the
  // G4SphericalSurface, then the solution must lie somewhere in between.

  G4int one_over_f = 8;  // one over fraction of a turn to go in each Step
  G4double turn_angle = 0.0;
  G4double dist_along = 0.0;
  G4double d_new;
  G4double fk = 1.0 / G4double( one_over_f );
  G4double scal = Scale();
  G4double d_old = HowNear( hx->position( dist_along ) );
  G4double rh = hx->GetRadius();    // radius of Helix
  G4Vector3D prp = hx->getPerp();   // perpendicular vector
  G4double prpmag = prp.mag();
  G4double rhp = rh / prpmag;
  G4int max_iter = one_over_f * HELIX_MAX_TURNS;

  // Take up to a user-settable number of turns along the Helix,
  // groping for an intersection point.

  for ( G4int k = 1; k < max_iter; k++ )
  {
    turn_angle = twopi * k / one_over_f;
    dist_along = turn_angle * std::fabs( rhp );
    d_new = HowNear( hx->position( dist_along ) );
    if ( ( d_old < 0.0 && d_new > 0.0 ) || ( d_old > 0.0 && d_new < 0.0 ) )
    {
      d_old = d_new;

      // Old and new points are on opposite sides of the G4SphericalSurface,
      // therefore a solution lies in between, use a binary search to pin the
      // point down to the surface precision, but don't do more than 50
      // iterations.

      G4int itr = 0;
      while ( std::fabs( d_new / scal ) > SURFACE_PRECISION )
      {
        itr++;
        if ( itr > 50 )  { return turn_angle; }
        turn_angle -= fk * pi;
        dist_along = turn_angle * std::fabs( rhp );
        d_new = HowNear( hx->position( dist_along ) );
        if (( d_old < 0.0 && d_new > 0.0 ) || ( d_old > 0.0 && d_new < 0.0 ))
          { fk *= -0.5; }
        else
          { fk *=  0.5; }
        d_old = d_new;
      }  //  end of while loop
      return turn_angle;  // this is the best solution
    }  //  end of if condition
  }  //  end of for loop

  //  Get here only if no solution is found, so return -1.0 to indicate that.

  return -1.0;
}
*/
