// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4CylindricalSurface.hh,v 1.4 2000-02-16 12:02:51 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
/*  /usr/local/gismo/repo/geometry/G4CylindricalSurface.h,v 1.16 1993/12/30 02:14:08 rensing Exp  */
//  File:  G4CylindricalSurface.h
//  Author:  Alan Breakstone

//  Contents ---------------------------------------------------------
//
//	G4CylindricalSurface
//
//  Description:
//   
//	C++ header file for Gismo G4CylindricalSurface class, derived from Surface class.
//	Uses the GmsListLink, G4ThreeVec, G4ThreeMat, Ray, Helix, and Surface
//	classes.
//
//  End --------------------------------------------------------------

//  Interface Dependencies -------------------------------------------

#ifndef __CYLINDER_H
#define __CYLINDER_H

#include "G4Surface.hh"

class G4ThreeMat;
 
//  End Interface Dependencies ---------------------------------------

//  Class  //

// class G4Surface;
class G4CylindricalSurface: public G4Surface
{
protected:          // make available to derived classes

  G4Vector3D axis;  // direction of axis of G4CylindricalSurface 
                    //   (unit vector)
  G4double radius;  // radius of G4CylindricalSurface

public:
  G4CylindricalSurface();
  G4CylindricalSurface( const G4Vector3D& o, 
			const G4Vector3D& a, 
			G4double r ); 
  
  virtual ~G4CylindricalSurface() {}
  
  //    G4CylindricalSurface( const G4CylindricalSurface& c ): 
  //                     G4Surface( c.origin )
  //			 { axis = c.axis;  radius = c.radius; }
  //
  
  G4String GetEntityType(){return G4String("Cylindrical_Surface");}
  
  virtual const char* NameOf() const { return "G4CylindricalSurface"; }
  
  virtual void PrintOn( G4std::ostream& os = G4cout ) const;
  
  int operator==( const G4CylindricalSurface& c )
  {
    return ( origin == c.origin    &&  
	     axis   == c.axis      && 
	     radius == c.radius       );
  }
  
  
  virtual G4double HowNear( const G4Vector3D& x ) const;
  
  //	virtual G4double distanceAlongRay( int which_way, const G4Ray* ry,
  //		G4Vector3D& p ) const;
  //	virtual G4double distanceAlongHelix( int which_way, 
  //                  const Helix* hx, G4Vector3D& p ) const;
  
  virtual G4Vector3D Normal( const G4Vector3D& p ) const;
  
  virtual G4Vector3D SurfaceNormal( const G4Point3D& p ) const;	
  
  virtual int Inside( const G4Vector3D& x ) const;
  
  virtual int WithinBoundary( const G4Vector3D& x ) const;
  
  virtual G4double Scale() const;

  //	virtual void rotate( G4double alpha, G4double beta, 
  //		G4double gamma, G4ThreeMat& m, int inverse ); 
  //	virtual void rotate( G4double alpha, G4double beta, 
  //		G4double gamma, int inverse ); 

  int Intersect(const G4Ray& ry);
  
  G4Vector3D GetAxis() const { return axis; }
  
  G4double GetRadius() const { return radius; }
  
  void SetRadius( G4double r );

private:

//	virtual G4double gropeAlongHelix( const Helix* hx ) const;
//
//
//  Description of functions -----------------------------------------
//
//  default constructor
//----->G4CylindricalSurface();
//
//  Normal constructor:first argument is the origin of the G4CylindricalSurface
//		       second argument is the axis of the G4CylindricalSurface
//		       third argument is the radius of the G4CylindricalSurface
//----->G4CylindricalSurface( const G4Vector3D& o, 
//                          const G4Vector3D& a, G4double r ); 
//
//  destructor 
//----->virtual ~G4CylindricalSurface() {}
//
//  copy constructor
//----->G4CylindricalSurface( const G4CylindricalSurface& c ): 
//                           Surface( c.origin )
//----->		 { axis = c.axis;  radius = c.radius; }
//
//  function to return class name
//----->virtual const char* NameOf() const { return "G4CylindricalSurface"; }
//
//  printing function
//----->virtual void PrintOn( G4std::ostream& os = G4cout ) const;
//
//  equality operator
//----->int operator==( const G4CylindricalSurface& c )
//----->	{ return origin == c.origin  &&  axis == c.axis
//----->	      && radius == c.radius; }
//
//  function which returns the distance from a point to a G4CylindricalSurface
//		the (input) argument is the point x
//		the distance is positive if the point is Inside,
//		negative if it is outside
//----->virtual G4double HowNear( const G4Vector3D& x ) const;
//
//  function which returns the distance along a Ray to enter or leave a
//  G4CylindricalSurface.  
//		the first (input) argument is +1 to leave or -1 to enter
//		the second (input) argument is a pointer to the Ray
//		the third (output) argument returns the intersection point
//----->virtual G4double distanceAlongRay( int which_way, const Ray* ry,
//----->	G4Vector3D& p ) const;
//
//  function which returns the distance along a Helix to enter or leave a
//  G4CylindricalSurface.  
//		the first (input) argument is +1 to leave or -1 to enter
//		the second (input) argument is a pointer to the Helix
//		the third (output) argument returns the intersection point
//----->virtual G4double distanceAlongHelix( int which_way, const Helix* hx,
//----->	G4Vector3D& p ) const;
//
//  function which returns the Normal unit vector to a 
//  G4CylindricalSurface at a point p
//  on (or nearly on) the G4CylindricalSurface
//----->virtual G4Vector3D Normal( const G4Vector3D& p ) const;
//
//  function which 
//          returns true (1) if the point x is Inside the G4CylindricalSurface,
//	    returns false (0) otherwise
//----->virtual int Inside( const G4Vector3D& x ) const;
//
//  function overwritten by finite-sized derived classes which returns
//		true (1) if the point x is within the boundary, false (0)
//		otherwise.
//		Since a G4CylindricalSurface is infinite in extent, the 
//              function will just check if the point is on the 
//              G4CylindricalSurface (to the surface precision).
//----->virtual int WithinBoundary( const G4Vector3D& x ) const;
//
//  function overwritten by finite-sized derived classes which returns
//		the radius, unless it is zero, in which case it returns
//		the smallest non-zero dimension. 
//		Used for Scale-invariant tests of surface thickness.
//----->virtual G4double Scale() const;
//
//  function to rotate the G4CylindricalSurface (4 input arguments) 
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
//  function to rotate the G4CylindricalSurface (4 input arguments) 
//		 first about global x-axis by angle alpha,
//		second about global y-axis by angle beta,
//		 third about global z-axis by angle gamma
//		the angles are assumed to be given in radians
//		the fourth (input) argument is an integer flag which if
//			non-zero reverses the order of the rotations
//----->virtual void rotate( G4double alpha, G4double beta, 
//----->	G4double gamma, int inverse ); 
//
//  functions to return the axis and radius of the G4CylindricalSurface
//----->direction GetAxis() const { return axis; }
//----->G4double GetRadius() const { return radius; }
//
//  function to change the radius of the G4CylindricalSurface
//----->void SetRadius( G4double r );
//
//
//  Private function to use a crude technique to find the intersection
//  of a Helix with a G4CylindricalSurface. It returns the turning angle along
//  the Helix at which the intersection occurs or -1.0 if no intersection
//  point is found.  The argument to the call is the pointer to the Helix.
//----->virtual G4double gropeAlongHelix( const Helix* hx ) const;
};

#endif


