// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ConicalSurface.hh,v 1.4 2000-02-16 12:02:51 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/*  /usr/local/gismo/repo/geometry/G4ConicalSurface.h,v 1.5 1993/12/30 02:13:59 rensing Exp  */
//  File:  G4ConicalSurface.h
//  Author:  Alan Breakstone

//  Contents ---------------------------------------------------------
//
//	G4ConicalSurface
//
//  Description:
//   
//	C++ header file for the Gismo G4ConicalSurface class, derived from 
//      Surface class.
//	Uses the GmsListLink, G4ThreeVec, G4ThreeMat, Ray, Helix, and Surface
//	classes.
//	A G4ConicalSurface is a semi-infinite conical surface defined by 
//      an axis and an opening angle, defined as the angle between the axis 
//      and the conical surface, with the origin being the apex of the cone.
//
//  End --------------------------------------------------------------

//  Interface Dependencies -------------------------------------------

#ifndef __CONICALSURFACE_H
#define __CONICALSURFACE_H

#include "G4Surface.hh"
class G4ThreeMat;
 
//  End Interface Dependencies ---------------------------------------

//  Class  //

class G4ConicalSurface: public G4Surface
{

private:
  G4Vector3D axis;   // direction of axis of G4ConicalSurface (unit vector)
  G4double   angle;  // half opening angle of G4ConicalSurface, in radians
		     // range is 0 < angle < PI/2

public:

  G4ConicalSurface();
  G4ConicalSurface( const G4Point3D& o, const G4Vector3D& a, G4double e ); 
  virtual ~G4ConicalSurface() {}

  G4String GetEntityType() { return G4String("Conical_Surface"); }

  // G4ConicalSurface( const G4ConicalSurface& c ): G4Surface( c.origin )
  //			 { axis = c.axis;  angle = c.angle; }
  
  virtual const char* NameOf() const { return "G4ConicalSurface"; }
  
  virtual void PrintOn( G4std::ostream& os = G4cout ) const;
  
  int operator==( const G4ConicalSurface& c )
  { 
    return origin == c.origin  &&  axis == c.axis && angle == c.angle; 
  }

  virtual G4double HowNear( const G4Vector3D& x ) const;

  //	virtual G4double distanceAlongRay( int which_way, const G4Ray* ry,
  //		G4Vector3D& p ) const;

  // Added 18.7-95
  void CalcBBox();
  
  // Added 18.7-95 , same as distanceAlongRay, but uses G4Ray.h
  int  Intersect( const G4Ray& ry );

  //	virtual G4double distanceAlongHelix( int which_way, 
  //                   const Helix* hx, G4Vector3D& p ) const;
  //    G4Vector3D Normal( const G4Vector3D& p ) const;

  virtual G4Vector3D SurfaceNormal( const G4Point3D& p ) const;	
  
  virtual int Inside( const G4Vector3D& x ) const;
  
  virtual int WithinBoundary( const G4Vector3D& x ) const;

  virtual G4double Scale() const { return 1.0; }

  //	virtual void rotate( G4double alpha, G4double beta, 
  //		G4double gamma, G4ThreeMat& m, int inverse ); 
  //	virtual void rotate( G4double alpha, G4double beta, 
  //		G4double gamma, int inverse ); 
  
  G4Vector3D GetAxis() const { return axis; }
  
  G4double GetAngle() const { return angle; }

  void SetAngle( G4double e );


private:

  //       	virtual G4double gropeAlongHelix( const Helix* hx ) const;
//
//  Description of functions -----------------------------------------
//
//  default constructor
//----->G4ConicalSurface();
//
//  Normal constructor:   first argument is the origin of the G4ConicalSurface
//			 second argument is the axis of the G4ConicalSurface
//			  third argument is the angle of the G4ConicalSurface
//----->G4ConicalSurface(const G4Vector3D& o, const G4Vector3D& a, G4double e);
//
//  destructor 
//----->virtual ~G4ConicalSurface() {}
//
//  copy constructor
//----->G4ConicalSurface( const G4ConicalSurface& c ): Surface( c.origin )
//----->		 { axis = c.axis;  angle = c.angle; }
//
//  function to return class name
//----->virtual const char* NameOf() const { return "G4ConicalSurface"; }
//
//  printing function
//----->virtual void PrintOn( G4std::ostream& os = G4cout ) const;
//
//  equality operator
//----->int operator==( const G4ConicalSurface& c )
//----->	{ return origin == c.origin  &&  axis == c.axis
//----->	      && angle == c.angle; }
//
//  function which returns the distance from a point to a G4ConicalSurface
//		the (input) argument is the point x
//		the distance is positive if the point is Inside,
//		negative if it is outside
//----->virtual G4double HowNear( const G4Vector3D& x ) const;
//
//  function which returns the distance along a Ray to enter or leave a
//  G4ConicalSurface.  
//		the first (input) argument is +1 to leave or -1 to enter
//		the second (input) argument is a pointer to the Ray
//		the third (output) argument returns the intersection point
//----->virtual G4double distanceAlongRay( int which_way, const Ray* ry,
//----->	G4Vector3D& p ) const;
//
//  function which returns the distance along a Helix to enter or leave a
//  G4ConicalSurface.  
//		the first (input) argument is +1 to leave or -1 to enter
//		the second (input) argument is a pointer to the Helix
//		the third (output) argument returns the intersection point
//----->virtual G4double distanceAlongHelix( int which_way, const Helix* hx,
//----->	G4Vector3D& p ) const;
//
//  function which returns the Normal unit vector to a G4ConicalSurface 
//  at a point p on (or nearly on) the G4ConicalSurface
//----->virtual G4Vector3D Normal( const G4Vector3D& p ) const;
//
//  function which returns 
//              true (1) if the point x is Inside the G4ConicalSurface,
//		false (0) otherwise
//----->virtual int Inside( const G4Vector3D& x ) const;
//
//  function overwritten by finite-sized derived classes which returns
//		true (1) if the point x is within the boundary, false (0)
//		otherwise.
//		Since a G4ConicalSurface is infinite in extent, the function 
//		will just check if the point is on the G4ConicalSurface 
//              (to the surface precision).
//----->virtual int WithinBoundary( const G4Vector3D& x ) const;
//
//  function overwritten by finite-sized derived classes which returns
//		a radius, unless it is zero, in which case it returns
//		the smallest non-zero dimension. 
//  		Since a semi-infinite cone has no Scale associated with it,
//		returns the arbitrary number 1.0.
//		Used for Scale-invariant tests of surface thickness.
//----->virtual G4double Scale() const { return 1.0; }
//
//  function to rotate the G4ConicalSurface (4 input arguments) 
//		 first about global x-axis by angle alpha,
//		second about global y-axis by angle beta,
//		 third about global z-axis by angle gamma
//		the angles are assumed to be given in radians
//		the fourth (output) argument gives the calculated rotation
//			matrix
//		the fifth (input) argument is an integer flag which if
//			non-zero reverses the order of the rotations
//----->virtual void rotate( G4double alpha, G4double beta, 
//----->	G4double gamma, G4ThreeMat& m, int inverse ); 
//
//  function to rotate the G4ConicalSurface (4 input arguments) 
//		 first about global x-axis by angle alpha,
//		second about global y-axis by angle beta,
//		 third about global z-axis by angle gamma
//		the angles are assumed to be given in radians
//		the fourth (input) argument is an integer flag which if
//			non-zero reverses the order of the rotations
//----->virtual void rotate( G4double alpha, G4double beta, 
//----->	G4double gamma, int inverse ); 
//
//  functions to return the axis and angle of the G4ConicalSurface
//----->direction GetAxis() const { return axis; }
//----->G4double GetAngle() const { return angle; }
//
//  function to change the angle of the G4ConicalSurface
//----->void SetAngle( G4double e );
//
//
//  Private function to use a crude technique to find the intersection
//  of a Helix with a G4ConicalSurface.  It returns the turning angle along the
//  Helix at which the intersection occurs or -1.0 if no intersection
//  point is found.  The argument to the call is the pointer to the Helix.
//----->virtual G4double gropeAlongHelix( const Helix* hx ) const;
};

#endif





