// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4CylindricalSurface.cc,v 1.3 2000-08-28 08:57:57 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4CylindricalSurface.cc
//
// ----------------------------------------------------------------------

#include "G4CylindricalSurface.hh"
#include "G4Sort.hh"
#include "G4Globals.hh"

G4CylindricalSurface::G4CylindricalSurface() : G4Surface()
{  
  // default constructor
  // default axis is ( 1.0, 0.0, 0.0 ), default radius is 1.0
  axis = G4Vector3D( 1.0, 0.0, 0.0 );
  radius = 1.0;
}


G4CylindricalSurface::G4CylindricalSurface( const G4Vector3D& o, 
					    const G4Vector3D& a,
					    G4double r ) //: G4Surface( o )
{ 
  // Normal constructor
  // require axis to be a unit vector
  G4double amag = a.mag();
  if ( amag != 0.0 )
    axis = a * (1/ amag);  // this makes the axis a unit vector
  else 
  {
    G4cerr << "Error in G4CylindricalSurface::G4CylindricalSurface--axis "
	   <<"has zero length\n"
	   << "\tDefault axis ( 1.0, 0.0, 0.0 ) is used.\n";
    
    axis = G4Vector3D( 1.0, 0.0, 0.0 );
  }

  //  Require radius to be non-negative
  if ( r >= 0.0 )
    radius = r;
  else 
  {
    G4cerr << "Error in G4CylindricalSurface::G4CylindricalSurface"
	   << "--asked for negative radius\n"
	   << "\tDefault radius of 1.0 is used.\n";
  
    radius = 1.0;
  }
	
  origin =o;
}


G4CylindricalSurface::~G4CylindricalSurface()
{
}

/*
G4CylindricalSurface::G4CylindricalSurface( const G4CylindricalSurface& c )
  : G4Surface( c.origin )
{
  axis = c.axis;  radius = c.radius;
}
*/

const char* G4CylindricalSurface::NameOf() const
{
  return "G4CylindricalSurface";
}

void G4CylindricalSurface::PrintOn( G4std::ostream& os ) const
{
  // printing function using C++ G4std::ostream class
  os << "G4CylindricalSurface surface with origin: " << origin << "\t"
     << "radius: " << radius << "\tand axis " << axis << "\n";
}

//G4int G4Surface::Intersect(const G4Ray& ry)
G4int G4CylindricalSurface::Intersect(const G4Ray& ry)
{

  //  L. Broglia : copy of G4FCylindricalSurface::Intersect
  
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

  // G4int which_way = -1; 
  //Originally a parameter.Read explanation above. 
  
  G4int which_way=1;
  
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
  G4int isoln   = 0, 
        maxsoln = 2;
  
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

  //  order the possible solutions by increasing distance along the Ray
  //  (G4Sorting routines are in support/G4Sort.h)
  sort_double( s, isoln, maxsoln-1 );

  //  now loop over each positive solution, keeping the first one (smallest
  //  distance along the Ray) which is within the boundary of the sub-shape
  //  and which also has the correct G4Vector3D with respect to the Normal to
  //  the G4CylindricalSurface at the intersection point
  for ( isoln = 0; isoln < maxsoln; isoln++ ) 
  {
    if ( s[isoln] >= kCarTolerance*0.5 ) 
    {
      if ( s[isoln] >= FLT_MAXX )  // quit if too large
	return 0;
      
      distance = s[isoln];
      closest_hit = ry.GetPoint( distance );
      G4double  tmp =  dhat * (Normal( closest_hit ));
      
      if ((tmp * which_way) >= 0.0 )
	if ( WithinBoundary( closest_hit ) == 1 )
	  distance =  distance*distance;
      
      return 1;
    }
  }
  
  //  get here only if there was no solution within the boundary, Reset
  //  distance and intersection point to large numbers
  distance = FLT_MAXX;
  closest_hit = lv;


  return 0;
}




G4double G4CylindricalSurface::HowNear( const G4Vector3D& x ) const
{
  //  Distance from the point x to the infinite G4CylindricalSurface.
  //  The distance will be positive if the point is Inside the 
  //  G4CylindricalSurface, negative if the point is outside.
  //  Note that this may not be correct for a bounded cylindrical object
  //  subclassed to G4CylindricalSurface.
  G4Vector3D d = x - origin;
  G4double dA = d * axis;
  G4double rad = sqrt( d.mag2() - dA*dA );
  G4double hownear = fabs( radius - rad );
 
  return hownear;
}


/*
G4double G4CylindricalSurface::distanceAlongRay( G4int which_way, const G4Ray* ry,
				   G4Vector3D& p ) const 
{  //  Distance along a Ray (straight line with G4Vector3D) to leave or enter
   //  a G4CylindricalSurface.  The input variable which_way should be set to +1
   //  to indicate leaving a G4CylindricalSurface, -1 to indicate entering a
   //  G4CylindricalSurface.
   //  p is the point of intersection of the Ray with the G4CylindricalSurface.
   //  If the G4Vector3D of the Ray is opposite to that of the Normal to
   //  the G4CylindricalSurface at the intersection point, it will not leave the
   //  G4CylindricalSurface.
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
	G4double Dist = FLT_MAXX;
	G4Vector3D lv ( FLT_MAXX, FLT_MAXX, FLT_MAXX );
	p = lv;
//  Origin and G4Vector3D unit vector of Ray.
	G4Vector3D x = ry->Position();
	G4Vector3D dhat = ry->Direction( 0.0 );
//  Axis unit vector of the G4CylindricalSurface.
	G4Vector3D ahat = GetAxis();
	G4int isoln = 0, maxsoln = 2;
//  array of solutions in distance along the Ray
//	G4double s[2] = { -1.0, -1.0 };
	G4double s[2];s[0] = -1.0; s[1]= -1.0 ;
//  calculate the two solutions (quadratic equation)
	G4Vector3D d = x - GetOrigin();
	G4double radius = GetRadius();
//  quit with no intersection if the radius of the G4CylindricalSurface is zero
	if ( radius <= 0.0 )
		return Dist;
	G4double dsq = d * d;
	G4double da = d * ahat;
	G4double dasq = da * da;
	G4double rsq = radius * radius;
	G4double qsq = dsq - dasq;
	G4double dira = dhat * ahat;
	G4double a = 1.0 - dira * dira;
	if ( a <= 0.0 )
		return Dist;
	G4double b = 2. * ( d * dhat - da * dira );
	G4double c = rsq - qsq;
	G4double radical = b * b + 4. * a * c; 
	if ( radical < 0.0 ) 
		return Dist;
	G4double root = sqrt( radical );
	s[0] = ( - b + root ) / ( 2. * a );
	s[1] = ( - b - root ) / ( 2. * a );
//  order the possible solutions by increasing distance along the Ray
//  (G4Sorting routines are in support/G4Sort.h)
	G4Sort_double( s, isoln, maxsoln-1 );
//  now loop over each positive solution, keeping the first one (smallest
//  distance along the Ray) which is within the boundary of the sub-shape
//  and which also has the correct G4Vector3D with respect to the Normal to
//  the G4CylindricalSurface at the intersection point
	for ( isoln = 0; isoln < maxsoln; isoln++ ) {
		if ( s[isoln] >= 0.0 ) {
			if ( s[isoln] >= FLT_MAXX )  // quit if too large
				return Dist;
			Dist = s[isoln];
			p = ry->Position( Dist );
			if ( ( ( dhat * Normal( p ) * which_way ) >= 0.0 )
			  && ( WithinBoundary( p ) == 1 ) )
				return Dist;
		}
	}
//  get here only if there was no solution within the boundary, Reset
//  distance and intersection point to large numbers
	p = lv;
	return FLT_MAXX;
}

*/


/*


G4double G4CylindricalSurface::distanceAlongHelix( G4int which_way,
                                   const Helix* hx,
				   G4Vector3D& p ) const 
{  //  Distance along a Helix to leave or enter a G4CylindricalSurface.
   //  The input variable which_way should be set to +1 to
   //  indicate leaving a G4CylindricalSurface, -1 to indicate entering a
   //  G4CylindricalSurface.
   //  p is the point of intersection of the Helix with the G4CylindricalSurface.
   //  If the G4Vector3D of the Helix is opposite to that of the Normal to
   //  the G4CylindricalSurface at the intersection point, it will not leave the
   //  G4CylindricalSurface.
   //  Similarly, if the G4Vector3D of the Helix is along that of the Normal 
   //  to the G4CylindricalSurface at the intersection point, it will not enter
   //  the G4CylindricalSurface.
   //  This method is called by all finite shapes sub-classed to
   //  G4CylindricalSurface.
   //  Use the virtual function table to check if the intersection point
   //  is within the boundary of the finite shape.
   //  If no valid intersection point is found, set the distance
   //  and intersection point to large numbers.
   //  Possible negative distance solutions are discarded.
	G4double Dist = FLT_MAXX;
	G4Vector3D lv ( FLT_MAXX, FLT_MAXX, FLT_MAXX );
	G4Vector3D zerovec;  // zero G4Vector3D
	p = lv;
	G4int isoln = 0, maxsoln = 4;
//  Array of solutions in turning angle
//	G4double s[4] = { -1.0, -1.0, -1.0, -1.0 };
	G4double s[4];s[0]=-1.0;s[1]= -1.0;s[2]= -1.0;s[3]= -1.0;


//  Flag set to 1 if exact solution is found
	G4int exact = 0;
//  Helix parameters
	G4double rh = hx->GetRadius();  	// radius of Helix
	G4Vector3D ah = hx->GetAxis();  	// axis of Helix
	G4Vector3D oh = hx->position();  // origin of Helix
	G4Vector3D dh = hx->direction( 0.0 ); // initial G4Vector3D of Helix
	G4Vector3D prp = hx->getPerp();	// perpendicular vector
	G4double prpmag = prp.Magnitude();
	G4double rhp = rh / prpmag;
//  G4CylindricalSurface parameters
	G4double rc = GetRadius();  // radius of G4CylindricalSurface
	if ( rc == 0.0 )   // quit if zero radius
		return Dist;
	G4Vector3D oc = GetOrigin();  // origin of G4CylindricalSurface
	G4Vector3D ac = GetAxis();  // axis of G4CylindricalSurface
//
//  Calculate quantities of use later on.
	G4Vector3D alpha = rhp * prp;
	G4Vector3D beta = rhp * dh;
	G4Vector3D gamma = oh - oc;
//  Declare variables used later on in several places.
	G4double rcd2 = 0.0, alpha2 = 0.0;
	G4double A = 0.0, B = 0.0, C = 0.0, F = 0.0, G = 0.0, H = 0.0;
	G4double CoverB = 0.0, radical = 0.0, root = 0.0, s1 = 0.0, s2 = 0.0;
	G4Vector3D ghat;
//
//  Set flag for special cases
	G4int special_case = 0;  // 0 means general case
//
//  Test to see if axes of Helix and G4CylindricalSurface are parallel, in which
//  case there are exact solutions.
	if ( ( fabs( ah.AngleBetween(ac) ) < FLT_EPSILO )
	  || ( fabs( ah.AngleBetween(ac) - M_PI ) < FLT_EPSILO ) ) {
	  	special_case = 1;
//  If, in addition, gamma is a zero vector or is parallel to the
//  G4CylindricalSurface axis, this simplifies the previous case.
		if ( gamma == zerovec ) {
			special_case = 3;
			ghat = gamma;
		}
		else {
			ghat = gamma / gamma.Magnitude();
			if ( ( fabs( ghat.AngleBetween(ac) ) < FLT_EPSILO )
			  || ( fabs( ghat.AngleBetween(ac) - M_PI ) <
			  				       FLT_EPSILO ) )
				special_case = 3;
		}
//  Test to see if, in addition to the axes of the Helix and G4CylindricalSurface
//  being parallel, the axis of the G4CylindricalSurface is perpendicular to the
//  initial G4Vector3D of the Helix.
		if ( fabs( ( ac * dh ) ) < FLT_EPSILO ) {
//  And, if, in addition to all this, the difference in origins of the Helix
//  and G4CylindricalSurface is  perpendicular to the initial G4Vector3D of the
//  Helix, there is a separate special case.
			if ( fabs( ( ghat * dh ) ) < FLT_EPSILO )
				special_case = 4;
		}
	}  // end of section with axes of Helix and G4CylindricalSurface parallel
//
//  Another peculiar case occurs if the axis of the G4CylindricalSurface and the
//  initial G4Vector3D of the Helix line up and their origins are the same.
//  This will require a higher order approximation than the general case.
	if ( ( ( fabs( dh.AngleBetween(ac) ) < FLT_EPSILO )
	    || ( fabs( dh.AngleBetween(ac) - M_PI ) < FLT_EPSILO ) )
	    && ( gamma == zerovec ) )
		special_case = 2;
//
//  Now all the special cases have been tagged, so solutions are found
//  for each case.  Exact solutions are indicated by setting exact = 1.
//  [For some reason switch doesn't work here, so do series of if's.]
	if ( special_case == 0  ) {  // approximate quadratic solutions
			A = beta * beta - ( beta * ac ) * ( beta * ac )
			  + gamma * alpha - ( gamma * ac ) * ( alpha * ac );
			B = 2.0 * gamma * beta 
			  - 2.0 * ( gamma * ac ) * ( beta * ac );
			C = gamma * gamma 
			  - ( gamma * ac ) * ( gamma * ac ) - rc * rc;
			if ( fabs( A ) < FLT_EPSILO ) {  // no quadratic term
				if ( B == 0.0 )  // no intersection, quit
					return Dist;
				else		 // B != 0
					s[0] = -C / B;
			}
			else {  // A != 0, general quadratic solution
				radical = B * B - 4.0 * A * C;
				if ( radical < 0.0 )  // no solution, quit 
					return Dist;
				root = sqrt( radical );
				s[0] = ( -B + root ) / ( 2.0 * A );
				s[1] = ( -B - root ) / ( 2.0 * A );
				if ( rh < 0.0 ) {
					s[0] = -s[0];
					s[1] = -s[1];
				}
				s[2] = s[0] + 2.0 * M_PI;
				s[3] = s[1] + 2.0 * M_PI;
			}
	}
//
	else if ( special_case == 1 ) {  // exact solutions
			exact = 1;
			H = 2.0 * ( alpha * alpha + gamma * alpha );
			F = gamma * gamma 
			  - ( ( gamma * ac ) * ( gamma * ac ) )
	             	  - rc * rc  + H;
			G = 2.0 * rhp * 
			    ( gamma * dh - ( gamma * ac ) * ( ac * dh ) );
			A = G * G  + H * H;
			B = -2.0 * F * H;
			C = F * F - G * G;
			if ( fabs( A ) < FLT_EPSILO ) {  // no quadratic term
				if ( B == 0.0 )  // no intersection, quit
					return Dist;
				else {  // B != 0
					CoverB = -C / B;
					if ( fabs( CoverB ) > 1.0 )
						return Dist;
					s[0] = acos( CoverB );
				}
			}
			else {  // A != 0, general quadratic solution
//  Try a different method of calculation using F, G, and H to avoid
//  precision problems.
//				radical = B * B - 4.0 * A * C;
//				if ( radical < 0.0 ) {
				  if ( fabs( H ) > FLT_EPSILO ) {
				    G4double r1 = G / H;
				    G4double r2 = F / H;
				    G4double radsq = 1.0 + r1*r1 - r2*r2;
				    if ( radsq < 0.0 )
				      return Dist;
				    root = G * sqrt( radsq );
				    G4double denominator = H * ( 1.0 + r1*r1 );
				    s1 = ( F + root ) / denominator;
				    s2 = ( F - root ) / denominator;
				  }
				  else
				    return Dist;
//				}  //  end radical < 0 condition
//				else {
//				  root = sqrt( radical );
//				  s1 = ( -B + root ) / ( 2.0 * A );
//				  s2 = ( -B - root ) / ( 2.0 * A );
//				}
				if ( fabs( s1 ) <= 1.0 ) {
					s[0] = acos( s1 );
					s[2] = 2.0 * M_PI - s[0];
				}
				if ( fabs( s2 ) <= 1.0 ) {
					s[1] = acos( s2 );
					s[3] = 2.0 * M_PI - s[1];
				}
//  Must take only solutions which satisfy original unsquared equation:
//  Gsin(s) - Hcos(s) + F = 0.  Take best solution of pair and set false
//  solutions to -1.  Only do this if the result is significantly different
//  from zero.
				G4double temp1 = 0.0, temp2 = 0.0;
				G4double rsign = 1.0;
				if ( rh < 0.0 ) rsign = -1.0;
				if ( s[0] > 0.0 ) {
				   temp1 =  G * rsign * sin( s[0] )
					  - H * cos( s[0] ) + F;
				   temp2 =  G * rsign * sin( s[2] )
					  - H * cos( s[2] ) + F;
				   if ( fabs( temp1 ) > fabs( temp2 ) )
					if ( fabs( temp1 ) > FLT_EPSILO ) 
					  s[0] = -1.0;
				   else
					if ( fabs( temp2 ) > FLT_EPSILO ) 
					  s[2] = -1.0;
				}
				if ( s[1] > 0.0 ) {
				   temp1 =  G * rsign * sin( s[1] )
					  - H * cos( s[1] ) + F;
				   temp2 =  G * rsign * sin( s[3] )
					  - H * cos( s[3] ) + F;
				   if ( fabs( temp1 ) > fabs( temp2 ) )
					if ( fabs( temp1 ) > FLT_EPSILO ) 
					  s[1] = -1.0;
				   else
					if ( fabs( temp2 ) > FLT_EPSILO ) 
					  s[3] = -1.0;
				}
			}
	}
//
	else if ( special_case == 2 ) {  // approximate solution
			G4Vector3D e = ah.cross( ac );
			G4double re = fabs( rhp ) * e.Magnitude();
			s[0] = sqrt( 2.0 * rc / re );
	}
//
	else if ( special_case == 3 ) {  // exact solutions
			exact = 1;
			alpha2 = alpha * alpha;
			rcd2 = rhp * rhp * ( 1.0 - ( (ac*dh) * (ac*dh) ) );
			A = alpha2 - rcd2;
			B = - 2.0 * alpha2;
			C = alpha2 + rcd2 - rc*rc;
			if ( fabs( A ) < FLT_EPSILO ) {  // no quadratic term
				if ( B == 0.0 )  // no intersection, quit
					return Dist;
				else {  // B != 0
					CoverB = -C / B;
					if ( fabs( CoverB ) > 1.0 )
						return Dist;
					s[0] = acos( CoverB );
				}
			}
			else {  // A != 0, general quadratic solution
				radical = B * B - 4.0 * A * C;
				if ( radical < 0.0 ) 
					return Dist;
				root = sqrt( radical );
				s1 = ( -B + root ) / ( 2.0 * A );
				s2 = ( -B - root ) / ( 2.0 * A );
				if ( fabs( s1 ) <= 1.0 )
					s[0] = acos( s1 );
				if ( fabs( s2 ) <= 1.0 )
					s[1] = acos( s2 );
			}
	}
//
	else if ( special_case == 4 ) {  // exact solution
			exact = 1;
			F = gamma * gamma 
			  - ( ( gamma * ac ) * ( gamma * ac ) )
	             	  - rc * rc;
			G = 2.0 * ( rhp * rhp  +  gamma * alpha );
			if ( G == 0.0 )  // no intersection, quit
				return Dist;
			G4double cs = 1.0 + ( F / G );
			if ( fabs( cs ) > 1.0 )  // no intersection, quit
				return Dist;
			s[0] = acos( cs );
	}
//
	else				// shouldn't get here
			return Dist;
//
// **************************************************************************
//
//  Order the possible solutions by increasing turning angle
//  (G4Sorting routines are in support/G4Sort.h).
	G4Sort_double( s, isoln, maxsoln-1 );
//
//  Now loop over each positive solution, keeping the first one (smallest
//  distance along the Helix) which is within the boundary of the sub-shape.
	for ( isoln = 0; isoln < maxsoln; isoln++ ) {
		if ( s[isoln] >= 0.0 ) {
//  Calculate distance along Helix and position and G4Vector3D vectors.
			Dist = s[isoln] * fabs( rhp );
			p = hx->position( Dist );
			G4Vector3D d = hx->direction( Dist );
			if ( exact == 0 ) {  //  only for approximate solns
//  Now do approximation to get remaining distance to correct this solution
//  iterate it until the accuracy is below the user-set surface precision.
			   G4double delta = 0.0;  
			   G4double delta0 = FLT_MAXX;
			   G4int dummy = 1;
			   G4int iter = 0;
			   G4int in0 = Inside( hx->position ( 0.0 ) );
			   G4int in1 = Inside( p );
			   G4double sc = Scale();
			   while ( dummy ) {
				iter++;
//  Terminate loop after 50 iterations and Reset distance to large number,
//  indicating no intersection with G4CylindricalSurface.
//  This generally occurs if the Helix curls too tightly to Intersect it.
				if ( iter > 50 ) {
					Dist = FLT_MAXX;
					p = lv;
					break;
				}
//  Find distance from the current point along the above-calculated
//  G4Vector3D using a Ray.
//  The G4Vector3D of the Ray and the Sign of the distance are determined
//  by whether the starting point of the Helix is Inside or outside of
//  the G4CylindricalSurface.
			   	in1 = Inside( p );
			   	if ( in1 ) {  //  current point Inside
				   if ( in0 ) {  //  starting point Inside
			   	      Ray* r = new Ray( p, d );
				      delta = 
				          distanceAlongRay( 1, r, p );
				      delete r;
				   }
				   else {       //  starting point outside
			   	      Ray* r = new Ray( p, -d );
				      delta = 
				         -distanceAlongRay( 1, r, p );
				      delete r;
				   }
				}
				else {        //  current point outside
				   if ( in0 ) {  //  starting point Inside
			   	      Ray* r = new Ray( p, -d );
				      delta = 
				         -distanceAlongRay( -1, r, p );
				      delete r;
				   }
				   else {        //  starting point outside
			   	      Ray* r = new Ray( p, d );
				      delta = 
				          distanceAlongRay( -1, r, p );
				      delete r;
				   }
				}
//  Test if distance is less than the surface precision, if so Terminate loop.
				if ( fabs( delta / sc ) <= SURFACE_PRECISION )
					break;
//  If delta has not changed sufficiently from the previous iteration, 
//  skip out of this loop.
				if ( fabs( ( delta - delta0 ) / sc ) <=
				                           SURFACE_PRECISION )
					break;
//  If delta has increased in absolute value from the previous iteration
//  either the Helix doesn't Intersect the G4CylindricalSurface or the approximate
//  solution is too far from the real solution.  Try groping for a solution.
//  If not found, Reset distance to large number, indicating no intersection with
//  the G4CylindricalSurface.
				if ( fabs( delta ) > fabs( delta0 ) ) {
					Dist = fabs( rhp ) * 
					       gropeAlongHelix( hx );
					if ( Dist < 0.0 ) {
						Dist = FLT_MAXX;
						p = lv;
					}
					else
						p = hx->position( Dist );
					break;
				}
//  Set old delta to new one.
				delta0 = delta;
//  Add distance to G4CylindricalSurface to distance along Helix.
				Dist += delta;
//  Negative distance along Helix means Helix doesn't Intersect
//  G4CylindricalSurface.
//  Reset distance to large number, indicating no intersection with
//  G4CylindricalSurface.
				if ( Dist < 0.0 ) {
					Dist = FLT_MAXX;
					p = lv;
					break;
				}
//  Recalculate point along Helix and the G4Vector3D.
				p = hx->position( Dist );
				d = hx->direction( Dist );
			   }  //  end of while loop
			}  //  end of exact == 0 condition
//  Now have best value of distance along Helix and position for this
//  solution, so test if it is within the boundary of the sub-shape
//  and require that it point in the correct G4Vector3D with respect to
//  the Normal to the G4CylindricalSurface.
			if ( ( Dist < FLT_MAXX ) &&
			     ( ( hx->direction( Dist ) * Normal( p ) *
			         which_way ) >= 0.0 ) &&
			     ( WithinBoundary( p ) == 1 ) )
				return Dist;
		}  // end of if s[isoln] >= 0.0 condition
	}  // end of for loop over solutions
//  if one gets here, there is no solution, so set distance along Helix
//  and position to large numbers
	Dist = FLT_MAXX;
	p = lv;
	return Dist;
}
*/


G4Vector3D G4CylindricalSurface::Normal( const G4Vector3D& p ) const
{ 
  //  return the Normal unit vector to the G4CylindricalSurface 
  //  at a point p on (or nearly on) the G4CylindricalSurface
 
  G4Vector3D n = ( p - origin ) - ( ( p - origin ) * axis ) * axis;
  G4double nmag = n.mag();
 
  if ( nmag != 0.0 )
    n = n * (1/nmag);
  
  return n;
}


G4Vector3D G4CylindricalSurface::SurfaceNormal( const G4Point3D& p ) const
{
  //  return the Normal unit vector to the G4CylindricalSurface at a point 
  //  p on (or nearly on) the G4CylindricalSurface
  
  G4Vector3D n = ( p - origin ) - ( ( p - origin ) * axis ) * axis;
  G4double nmag = n.mag();
  
  if ( nmag != 0.0 )
    n = n * (1/nmag);
  
  return n;
}


G4int G4CylindricalSurface::Inside ( const G4Vector3D& x ) const
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


G4int G4CylindricalSurface::WithinBoundary( const G4Vector3D& x ) const
{ 
  //  return 1 if point x is on the G4CylindricalSurface, otherwise return zero
  //  base this on the surface precision factor set in support/globals.h
  if ( fabs( HowNear( x ) / Scale() ) <= SURFACE_PRECISION )
    return 1;
  else
    return 0;
}


G4double G4CylindricalSurface::Scale() const
{
  //  Returns the radius of a G4CylindricalSurface unless it is zero, in which
  //  case returns the arbitrary number 1.0.
  //  This is ok since derived finite-sized classes will overwrite this.
  //  Used for Scale-invariant tests of surface thickness.
  if ( radius == 0.0 )
    return 1.0;
  else
    return radius;
}


//void G4CylindricalSurface::rotate( G4double alpha, G4double beta,
//		       G4double gamma, G4ThreeMat& m, G4int inverse )
// //  rotate G4CylindricalSurface  first about global x-axis by angle alpha,
//                  second about global y-axis by angle beta,
//               and third about global z-axis by angle gamma
//  by creating and using G4ThreeMat objects in Surface::rotate
//  angles are assumed to be given in radians
//  if inverse is non-zero, the order of rotations is reversed
//  the axis is rotated here, the origin is rotated by calling
//  Surface::rotate
//	G4Surface::rotate( alpha, beta, gamma, m, inverse );
//	axis = m * axis;
//}

//void G4CylindricalSurface::rotate( G4double alpha, G4double beta,
//		       G4double gamma, G4int inverse )
//{  //  rotate G4CylindricalSurface  first about global x-axis by angle alpha,
//                  second about global y-axis by angle beta,
//               and third about global z-axis by angle gamma
//  by creating and using G4ThreeMat objects in Surface::rotate
//  angles are assumed to be given in radians
//  if inverse is non-zero, the order of rotations is reversed
//  the axis is rotated here, the origin is rotated by calling
//  Surface::rotate
//	G4ThreeMat m;
//	G4Surface::rotate( alpha, beta, gamma, m, inverse );
//	axis = m * axis;
//}


void G4CylindricalSurface::SetRadius( G4double r )
{ 
  //  Reset the radius of the G4CylindricalSurface
  //  Require radius to be non-negative
  if ( r >= 0.0 )
    radius = r;
  //  use old value (do not change radius) if out of the range, 
  //  but Print message
  else 
  {
    G4cerr << "Error in G4CylindricalSurface::SetRadius"
	   << "--asked for negative radius\n"
	   << "\tDefault radius of " << radius << " is used.\n";
  }
}


/*
G4double G4CylindricalSurface::gropeAlongHelix( const Helix* hx ) const
{  //  Grope for a solution of a Helix intersecting a G4CylindricalSurface.
   //  This function returns the turning angle (in radians) where the
   //  intersection occurs with only positive values allowed, or -1.0 if
   //  no intersection is found.
//  The idea is to start at the beginning of the Helix, then take steps
//  of some fraction of a turn.  If at the end of a Step, the current position
//  along the Helix and the previous position are on opposite sides of the
//  G4CylindricalSurface, then the solution must lie somewhere in between.
	G4int one_over_f = 8;  // one over fraction of a turn to go in each Step
	G4double turn_angle = 0.0;
	G4double dist_along = 0.0;
	G4double d_new;
	G4double fk = 1.0 / G4double( one_over_f );
	G4double scal = Scale();
	G4double d_old = HowNear( hx->position( dist_along ) );
	G4double rh = hx->GetRadius();	// radius of Helix
	G4Vector3D prp = hx->getPerp();	// perpendicular vector
	G4double prpmag = prp.Magnitude();
	G4double rhp = rh / prpmag;
	G4int max_iter = one_over_f * HELIX_MAX_TURNS;
//  Take up to a user-settable number of turns along the Helix,
//  groping for an intersection point.
	for ( G4int k = 1; k < max_iter; k++ ) {
		turn_angle = 2.0 * M_PI * k / one_over_f;
		dist_along = turn_angle * fabs( rhp );
		d_new = HowNear( hx->position( dist_along ) );
		if ( ( d_old < 0.0 && d_new > 0.0 ) ||
		     ( d_old > 0.0 && d_new < 0.0 ) ) {
			d_old = d_new;
//  Old and new points are on opposite sides of the G4CylindricalSurface, therefore
//  a solution lies in between, use a binary search to pin the point down
//  to the surface precision, but don't do more than 50 iterations.
			G4int itr = 0;
			while ( fabs( d_new / scal ) > SURFACE_PRECISION ) {
				itr++;
				if ( itr > 50 )
					return turn_angle;
				turn_angle -= fk * M_PI;
				dist_along = turn_angle * fabs( rhp );
				d_new = HowNear( hx->position( dist_along ) );
				if ( ( d_old < 0.0 && d_new > 0.0 ) ||
				     ( d_old > 0.0 && d_new < 0.0 ) )
					fk *= -0.5;
				else
					fk *=  0.5;
				d_old = d_new;
			}  //  end of while loop
			return turn_angle;  // this is the best solution
		}  //  end of if condition
	}  //  end of for loop
//  Get here only if no solution is found, so return -1.0 to indicate that.
	return -1.0;
}
*/
