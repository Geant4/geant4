// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FCylindricalSurface.hh,v 1.8 2000-08-28 08:57:45 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4FCylindricalSurface
//
// Class Description:
//   
// Definition of a generic bounded cylindrical surface.
//
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

// The code for G4CylindricalSurface has been derived from the original
// implementation in the "Gismo" package.
//
// Author: A.Breakstone
// Adaptation: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef __FCYLINDER_H
#define __FCYLINDER_H

#include "G4PointRat.hh"
#include "G4Axis2Placement3D.hh"
#include "G4Surface.hh"


class G4FCylindricalSurface : public G4Surface
{

 public:  // with description

  G4FCylindricalSurface();
    // Default constructor.

  G4FCylindricalSurface( const G4Point3D&  o, 
			 const G4Vector3D& a,
			 G4double    r, 
			 G4double    l );
    // o : origin.
    // a : axis.
    // r : radius.
    // l : length.

  virtual ~G4FCylindricalSurface();
    // Destructor.

  G4FCylindricalSurface(const G4FCylindricalSurface& c);
    // Copy constructor.

  inline G4int operator==( const G4FCylindricalSurface& c );
    // Equality operator.


  virtual G4Vector3D SurfaceNormal( const G4Point3D& p ) const;
    // Returns the Normal unit vector to the G4FCylindricalSurface at a point 
    // p on (or nearly on) the G4FCylindricalSurface.

  virtual G4int Inside( const G4Vector3D& x ) const;
    // Returns 1 if the point x is Inside the G4FCylindricalSurface,
    // returns 0 otherwise.

  inline G4String GetEntityType();
    // Returns the shape identifier.

  G4int Intersect(const G4Ray&);	
    // Counts the number of intersections of a bounded cylindrical surface
    // by a ray. At first, it calculates the intersections with the infinite 
    // cylindrical surface. Then, it counts the intersections within the
    // finite cylindrical surface boundaries, and sets the "distance" to the 
    // closest distance from the start point to the nearest intersection.
    // If the point is on the surface it returns either the intersection with
    // the opposite surface or kInfinity.
    // If no intersection is found, it sets distance = kInfinity and returns 0.

  virtual G4double HowNear( const G4Vector3D& x ) const;
    // Computes the shortest distance from the point x to the
    // G4FCylindricalSurface. The distance will be always positive. 

  void CalcBBox();
    // Computes the bounding box.

  virtual const char* NameOf() const; 
    // Returns the class name.

  virtual void PrintOn( G4std::ostream& os = G4cout ) const;
    // Printing function.

  virtual G4int WithinBoundary( const G4Vector3D& x ) const;
    // Returns 1 if the point x is within the boundary, 0 otherwise.

  virtual G4double Scale() const;
    // Return the radius of a G4FCylindricalSurface.
    // Used for Scale-invariant tests of surface thickness.
    // If the radius is zero, it returns the length.

  virtual G4double Area() const;
    // Calculates the area of a G4FCylindricalSurface.
  
  virtual void resize( G4double r, G4double l );
    // Changes the radius and length of the G4FCylindricalSurface:
    // the first (input) argument is the new radius;
    // the second (input) argument is the new length.

  inline G4double GetLength() const;
  inline G4Vector3D GetAxis() const;
  inline G4double GetRadius() const;
  void SetRadius( G4double r );
    // Get/Set methods for the geometrical data of the G4FCylindricalSurface.

 public:  // without description

  void InitValues();
    // Re-calculates the private values of the G4FCylindrical surface
    // before the Intersect and HowNear function if the G4FCylindrical
    // was created by the STEP interface.

 protected:

  G4Axis2Placement3D Position;
  G4double  radius; 
  G4double  length;

};

#include "G4FCylindricalSurface.icc"

#endif
