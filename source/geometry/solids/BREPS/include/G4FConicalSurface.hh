// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FConicalSurface.hh,v 1.6 2000-02-16 12:02:52 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef __FCONIC_H
#define __FCONIC_H

#include "G4PointRat.hh"
#include "G4Axis2Placement3D.hh"
#include "G4Surface.hh"


//     Position.axis|
//                  |
//    --         ---|---   small_radius  
//  l  |        /   |   \
//  e  |       /    |    \
//  n  |      /     |     \
//  g  |     /      |      \
//  t  |    /       |       \
//  h  |   /        |        \
//    --   ---------|---------  large_radius
//               Position


class G4FConicalSurface: public G4Surface //: public G4ConicalSurface
{
protected:

  G4double length;    	// length of G4FConicalSurface
  G4double small_radius;// small radius of G4FConicalSurface, can be zero
  G4double large_radius;// large radius of G4FConicalSurface, must be 
                        // greater than the small radius
                        // Note that the angle of the G4ConicalSurface is
                        // calculated from these three quantities.
	
  G4Axis2Placement3D Position;

  // Add by L. Broglia
  G4double tan_angle;

public:
  
  G4FConicalSurface() //: G4ConicalSurface() 
  {
    length       = 1.0;
    small_radius = 0.0;
    large_radius = 1.0;
    
    // Add by L. Broglia
    tan_angle = (large_radius-small_radius)/length;
  }
	
// constructor utilized into G4BREPSolidPCone
//  default constructor
//----->G4FConicalSurface() : G4ConicalSurface() { length = 1.0;
//----->		     small_radius = 0.0;
//----->		     large_radius = 1.0; }
//
//  Normal constructor:  first  argument is the origin of the G4FConicalSurface
//			 second argument is the axis of the G4FConicalSurface
//		         third  argument is the length of the G4FConicalSurface
//		         fourth argument is the small radius of the 
//                                          G4FConicalSurface
//			 fifth argument is the large radius of the 
//                                          G4FConicalSurface

  G4FConicalSurface( const G4Point3D& o, const G4Vector3D& a,
		     G4double l, G4double sr, G4double lr );

  G4FConicalSurface( const G4FConicalSurface& c );
  
  ~G4FConicalSurface() {}

  virtual G4Vector3D SurfaceNormal( const G4Point3D& p ) const;	

// Return 0 if point x is outside G4ConicalSurface, 1 if Inside.
  virtual int Inside( const G4Vector3D& x ) const;

  G4String GetEntityType(){return G4String("FConical_Surface");}

  // STEP additions
//
//  function to return class name
  virtual const char* Name() const { return "G4FConicalSurface"; }
//  printing function
  virtual void PrintOn( G4std::ostream& os = G4cout ) const;

//  equality operator
  int operator==( const G4FConicalSurface& c );

// This function count the number of intersections of a 
// bounded conical surface by a ray.
// At first, calculates the intersections with the semi-infinite 
// conical surfsace. After, count the intersections within the
// finite conical surface boundaries, and set "distance" to the 
// closest distance from the start point to the nearest intersection
// If the point is on the surface it returns or the intersection with
// the opposite surface or kInfinity
// If no intersection is founded, set distance = kInfinity and
// return 0

  int  Intersect( const G4Ray& ry ) ;
  void CalcBBox();

  // Add by L. Broglia

// Shortest distance from the point x to the G4FConicalSurface.
// The distance will be positive always positive

 virtual G4double HowNear( const G4Vector3D& x ) const;
  
//  function which returns true (1) if the point x is within the boundary
//		returns false (0) otherwise

  virtual int WithinBoundary( const G4Vector3D& x ) const;

//  function to return the size of a G4FConicalSurface.
//		Used for Scale-invariant tests of surface thickness.
//		If the small radius is zero, returns the large radius.

  virtual G4double Scale() const;

//  function to calculate the Area of a G4FConicalSurface
  virtual G4double Area() const;

//  function to change the radii and length of the G4FConicalSurface
//		the first (input) argument is the new length
//		the second (input) argument is the new small radius
//		the third (input) argument is the new large radius

  virtual void resize( G4double l, G4double sr, G4double lr );
  
//  functions to return the dimensions of the G4FConicalSurface

  G4double GetLength()      const { return length;       }

//  functions to return the dimensions of the G4FConicalSurface

  G4double GetSmallRadius() const { return small_radius; }

//  functions to return the dimensions of the G4FConicalSurface

  G4double GetLargeRadius() const { return large_radius; }

//  functions to return the dimensions of the G4FConicalSurface

  G4double GetTan_Angle()   const { return tan_angle; }

};

#endif









