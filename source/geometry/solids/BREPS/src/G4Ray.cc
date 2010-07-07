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
// $Id: G4Ray.cc,v 1.13 2010-07-07 14:45:31 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4Ray.cc
//
// ----------------------------------------------------------------------

#include "G4Ray.hh"
#include "G4PointRat.hh"

G4Ray::G4Ray()
  : r_min(0.), r_max(0.)
{
}

G4Ray::G4Ray(const G4Point3D& start0, const G4Vector3D& dir0)
{
  Init(start0, dir0);
}

G4Ray::~G4Ray()
{
}


const G4Plane& G4Ray::GetPlane(G4int number_of_plane) const
{
  if(number_of_plane==1)
    { return plane2; }
  else
    { return plane1; }
}


void G4Ray::CreatePlanes()
{
  // Creates two orthogonal planes(plane1,plane2) the ray (rray) 
  // situated in the intersection of the planes. The planes are 
  // used to project the surface (nurb) in two dimensions.
  
  G4Vector3D RayDir    = dir;
  G4Point3D  RayOrigin = start;

  G4Point3D  p1, p2, p3, p4;
  G4Vector3D dir1, dir2;
  G4Vector3D invdir = G4Vector3D( PINFINITY );
  
  if(!NearZero(RayDir.x(), SQRT_SMALL_FASTF)) 
    { invdir.setX(1.0 / RayDir.x()); }
    
  if(!NearZero(RayDir.y(), SQRT_SMALL_FASTF)) 
    { invdir.setY(1.0 / RayDir.y()); }
    
  if(!NearZero(RayDir.z(), SQRT_SMALL_FASTF)) 
    { invdir.setZ(1.0 / RayDir.z()); }

  MatVecOrtho(dir1, RayDir);
  
  Vcross( dir2, RayDir, dir1);
  Vmove(p1, RayOrigin);
  Vadd2(p2, RayOrigin, RayDir);
  Vadd2(p3, RayOrigin, dir1);
  Vadd2(p4, RayOrigin, dir2);

  CalcPlane3Pts( plane1, p1, p3, p2);
  CalcPlane3Pts( plane2, p1, p2, p4);
}


void G4Ray::MatVecOrtho(register G4Vector3D &out,
                        register const G4Vector3D &in )
{
  register G4double f;
  G4int             i_Which;

  if( NearZero(in.x(), 0.0001)
   && NearZero(in.y(), 0.0001)
   && NearZero(in.z(), 0.0001) )  
  {
    Vsetall( out, 0 );
    return;
  }
 
  //      Find component closest to zero 
  f = std::fabs(in.x());
  i_Which=0;
  
  if( std::fabs(in.y()) < f )
  {
    f = std::fabs(in.y());
    i_Which=1;
  }
  
  if( std::fabs(in.z()) < f )
  {
    i_Which=2;
  }
  
  if(!i_Which)
  {
    f = std::sqrt((in.y())*(in.y())+(in.z())*(in.z()));    // hypot(in.y(),in.z())
  }
  else
  {
    if(i_Which==1)
    {
      f = std::sqrt((in.z())*(in.z())+(in.x())*(in.x()));  // hypot(in.z(),in.x())
    }
    else
    {
      f = std::sqrt((in.x())*(in.x())+(in.y())*(in.y()));  // hypot(in.x(),in.y())
    }
  }
  if( NearZero( f, SMALL ) )
  {
    Vsetall( out, 0 );
    return;
  }
    
  f = 1.0/f;
    
  if(!i_Which)
  {
    out.setX(0.0);
    out.setY(-in.z()*f);
    out.setZ( in.y()*f);
  }
  else
  {
    if(i_Which==1)
    {
      out.setY(0.0);
      out.setZ(-in.x()*f);
      out.setX( in.y()*f);
    }
    else
    {
      out.setZ(0.0);
      out.setX(-in.z()*f);
      out.setY( in.y()*f);
    }
  }
} 


//    			CALC_PLANE_3PTS
//
//  Find the equation of a G4Plane that contains three points.
//  Note that Normal vector created is expected to point out (see vmath.h),
//  so the vector from A to C had better be counter-clockwise
//  (about the point A) from the vector from A to B.
//  This follows the outward-pointing Normal convention, and the
//  right-hand rule for cross products.
//
/*
                        C
                        *
                        |\
                        | \
           ^     N      |  \
           |      \     |   \
           |       \    |    \
           |C-A     \   |     \
           |         \  |      \
           |          \ |       \
                       \|        \
                        *---------*
                        A         B
                           ----->
                            B-A
*/
//  If the points are given in the order A B C (eg, *counter*-clockwise),
//  then the outward pointing surface Normal N = (B-A) x (C-A).
//
//  Explicit Return -
//       0      OK
//      -1      Failure.  At least two of the points were not distinct,
//              or all three were colinear.
//
//  Implicit Return -
//      G4Plane   The G4Plane equation is stored here.


G4int G4Ray::CalcPlane3Pts(G4Plane &plane1,
			   const G4Point3D& a,
			   const G4Point3D& b,
			   const G4Point3D& c )
{
  // Creates the two orthogonal planes which are needed in projecting the
  // surface into 2D.

  G4Vector3D	B_A;
  G4Vector3D    C_A;
  G4Vector3D    C_B;
  
  register G4double mag;

  Vsub2( B_A, b, a );
  Vsub2( C_A, c, a );
  Vsub2( C_B, c, b );

  Vcross( plane1, B_A, C_A );

  //	Ensure unit length Normal 
  mag = Magnitude(plane1);
  if( mag  <= SQRT_SMALL_FASTF )
  {
    return(-1);//	 FAIL 
  }
  
  mag = 1/mag;
  
  G4Plane pl2(plane1);
  Vscale( plane1, pl2, mag );

  //     Find distance from the origin to the G4Plane 
  plane1.d = Vdot( plane1, a );
  
  return(0);	//ok
}


void G4Ray::RayCheck()
{
  // Check that the ray has a G4Vector3D...
  if (dir==G4Vector3D(0, 0, 0)) 
  {
    G4Exception("G4Ray::RayCheck()", "InvalidInput", FatalException,
                "Invalid zero direction given !");
  }

  // Make sure that the vector is unit length
  dir= dir.unit();
  r_min = 0;
  r_max = 0;
}


void G4Ray::Vcross(G4Plane &a, const G4Vector3D &b, const G4Vector3D &c) 
{ 
  a.a = b.y()  * c.z()  - b.z()  * c.y() ;
  a.b = b.z()  * c.x()  - b.x()  * c.z() ;
  a.c = b.x()  * c.y()  - b.y()  * c.x() ;
}


void G4Ray::Vcross(G4Vector3D &a, const G4Vector3D &b, const G4Vector3D &c) 
{ 
  a.setX(b.y()  * c.z()  - b.z()  * c.y()) ;
  a.setY(b.z()  * c.x()  - b.x()  * c.z()) ;
  a.setZ(b.x()  * c.y()  - b.y()  * c.x()) ;
}


void G4Ray::Vmove(G4Point3D &a, const G4Point3D &b) 
{ 
  a.setX(b.x());
  a.setY(b.y());
  a.setZ(b.z());
}


void G4Ray::Vadd2(G4Point3D &a, const G4Point3D &b, const G4Vector3D &c) 
{
  a.setX(b.x() + c.x()) ;
  a.setY(b.y() + c.y()) ;
  a.setZ(b.z() + c.z()) ;
}	
  

void G4Ray::Vsub2(G4Vector3D &a, const G4Point3D &b, const G4Point3D &c) 
{
  a.setX(b.x() - c.x());
  a.setY(b.y() - c.y());
  a.setZ(b.z() - c.z());
}


void G4Ray::Vscale(G4Plane& a, const G4Plane& b, G4double c) 
{ 
  a.a = b.a * c;
  a.b = b.b * c;
  a.c = b.c * c;
}


G4double G4Ray::Vdot(const G4Plane &a, const G4Point3D &b) 
{
  return (a.a * b.x() + 
	  a.b * b.y() + 
	  a.c * b.z());
}
  

G4double G4Ray::Magsq(const G4Plane &a) 
{
  return ( a.a * a.a + a.b * a.b + a.c *a.c );
}
  

G4double G4Ray::Magnitude(const G4Plane &a) 
{
  return (std::sqrt( Magsq( a )) );
}
