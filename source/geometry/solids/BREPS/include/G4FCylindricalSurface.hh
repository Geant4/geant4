// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FCylindricalSurface.hh,v 1.4 1999-05-19 16:55:46 magni Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef __FCYLINDER_H
#define __FCYLINDER_H

#include "G4PointRat.hh"
#include "G4Axis2Placement3D.hh"
#include "G4Surface.hh"


//     Position.axis|        radius
//                 >|---|<---------
//                  |   
//    --        +---|---+  
//  l  |        |   |   |
//  e  |        |   |   | 
//  n  |        |   |   |  
//  g  |        |   |   |   
//  t  |        |   |   |    
//  h  |        |   |   |     
//    --        +---|---+  
//               Position
 


class G4FCylindricalSurface: public G4Surface
{
 protected:

  G4Axis2Placement3D Position;
  G4double  radius; 
  G4double  length;


 public:

  //  default constructor
  G4FCylindricalSurface()
  {
    length = 1.0; 
  }
  
  // Constructor utilized by the G4BREPSolidPCone
  //  first argument  is the origin 
  //  second argument is the axis  
  //  third argument  is the radius 
  //  fourth argument is the length 
  G4FCylindricalSurface( const G4Point3D&  o, 
			 const G4Vector3D& a,
			 const G4double    r, 
			 const G4double    l );

  //  destructor 
  ~G4FCylindricalSurface() {}
  
  //  copy constructor
  G4FCylindricalSurface(const G4FCylindricalSurface& c);

  virtual G4Vector3D SurfaceNormal( const G4Point3D& p ) const;
  
  virtual int Inside( const G4Vector3D& x ) const;
 
  // 
  G4String GetEntityType() 
  {
    return G4String("Cylindrical_Surface");
  }
  
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
  int Intersect(const G4Ray&);	

  // Shortest distance from the point x to the G4FCylindricalSurface.
  // The distance will be always positive. 
   virtual G4double HowNear( const G4Vector3D& x ) const;

  //
  void CalcBBox();
  
  //  function to return class name   
  virtual char *NameOf() const 
  {
    return "G4FCylindricalSurface"; 
  }
  
  //  printing function
  virtual void PrintOn( ostream& os = G4cout ) const;
  
  //  equality operator
  int operator==( const G4FCylindricalSurface& c );
  
  //  function which returns true  (1) if the point x is within the boundary
  //	             returns false (0) otherwise
  virtual int WithinBoundary( const G4Vector3D& x ) const;
  
  //  function to return the radius of a G4FCylindricalSurface.
  //	 Used for Scale-invariant tests of surface thickness.
  //	 If the radius is zero, returns the length.
  virtual G4double Scale() const;
  
  //  function to calculate the Area of a G4FCylindricalSurface
  virtual G4double Area() const 
  {
    return ( 2.0 * M_PI * radius * length );
  }
  
  //  function to change the radius and length of the G4FCylindricalSurface
  //	  the first (input) argument is the new radius
  //	  the second (input) argument is the new length
  virtual void resize( G4double r, G4double l );
  
  //  function to return the length of the G4FCylindricalSurface
  G4double GetLength() const 
  { 
    return length; 
  }
    
  G4Vector3D GetAxis() const { return Position.GetAxis(); }
  
  G4double GetRadius() const { return radius; }
  
  void SetRadius( G4double r );

  // Re-calculates the private values of the G4FCylindrical surface
  // before the Intersect and HowNear function if the G4FCylindrical
  // was created by the STEP interface
  void InitValues();
};

#endif






