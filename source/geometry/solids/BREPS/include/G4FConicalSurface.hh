// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FConicalSurface.hh,v 1.1 1999-01-07 16:07:31 gunter Exp $
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
			     
  G4FConicalSurface( const G4Point3D& o, const G4Vector3D& a,
		     G4double l, G4double sr, G4double lr );

  G4FConicalSurface( const G4FConicalSurface& c );
  
  ~G4FConicalSurface() {}

  virtual G4Vector3D SurfaceNormal( const G4Point3D& p ) const;	

  virtual int Inside( const G4Vector3D& x ) const;

  G4String GetEntityType(){return G4String("FConical_Surface");}

  // STEP additions
  virtual char *Name() const { return "G4FConicalSurface"; }
  virtual void PrintOn( ostream& os = G4cout ) const;

  int operator==( const G4FConicalSurface& c );

  int  Intersect( const G4Ray& ry ) ;
  void CalcBBox();

  // Add by L. Broglia
  virtual G4double HowNear( const G4Vector3D& x ) const;
  
  inline void Comp( G4Vector3D& v, G4Point3D& min , G4Point3D& max)
  {
    if(v.x() > max.x() ) max.setX(v.x());
    if(v.y() > max.y() ) max.setY(v.y());
    if(v.z() > max.z() ) max.setZ(v.z());

    if(v.x() < min.x()) min.setX(v.x());
    if(v.y() < min.y()) min.setY(v.y());
    if(v.z() < min.z()) min.setZ(v.z()); 
  }

  
  virtual int WithinBoundary( const G4Vector3D& x ) const;
  virtual G4double Scale() const;
  virtual G4double Area() const;
  virtual void resize( G4double l, G4double sr, G4double lr );
  
  G4double GetLength()      const { return length;       }
  G4double GetSmallRadius() const { return small_radius; }
  G4double GetLargeRadius() const { return large_radius; }
  G4double GetTan_Angle()   const { return tan_angle; }

//  Description of functions -----------------------------------------
//
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
//----->G4FConicalSurface( const G4ThreeVec& o, const G4ThreeVec& a,
//----->	        G4double l, G4double sr, G4double lr );
//
//  destructor 
//----->virtual ~G4FConicalSurface() {}
//
//  copy constructor
//----->G4FConicalSurface( const G4FConicalSurface& c );
//
//  function to return class name
//----->virtual char *NameOf() const { return "G4FConicalSurface"; }
//
//  printing function
//----->virtual void PrintOn( ostream& os = G4cout ) const;
//
//  equality operator
//----->int operator==( const G4FConicalSurface& c );
//
//  function which returns true (1) if the point x is within the boundary
//		returns false (0) otherwise
//----->virtual int WithinBoundary( const G4ThreeVec& x ) const;
//
//  function to return the size of a G4FConicalSurface.
//		Used for Scale-invariant tests of surface thickness.
//		If the small radius is zero, returns the large radius.
//----->virtual G4double Scale() const;
//
//  function to calculate the Area of a G4FConicalSurface
//----->virtual G4double Area() const; 
//
//  function to change the radii and length of the G4FConicalSurface
//		the first (input) argument is the new length
//		the second (input) argument is the new small radius
//		the third (input) argument is the new large radius
//----->virtual void resize( G4double l, G4double sr, G4double lr );
//
//  functions to return the dimensions of the G4FConicalSurface
//----->G4double GetLength() const { return length; }
//----->G4double GetSmallRadius() const { return small_radius; }
//----->G4double GetLargeRadius() const { return large_radius; }
};

#endif









